prep_Kaplan_Meier_data <- function (step, taskname) {
  
  ll_bigF <- if (taskname == "AON") AON_bigF else leaky_bigF
  
  # read in data (=observations)
  if (step == "A") {
    dt_observed <- 
      read_data_A(step) %>% 
      transmute(
        Location,
        Drug_regimen,
        Day,
        m  = cum_risk,
        q1 = cum_risk - sd_risk,
        q2 = cum_risk + sd_risk
      ) %>%
      mutate(Data = "Clinical trial observations")
  } else if (step == "B") {
    dt_observed <-
      read_csv(file.path(dirX_inp(step), "data.csv")) %>%
      select(Drug_regimen, Day, m, q1, q2) %>%
      mutate(
        m  = 100 * m,
        q1 = 100 * q1,
        q2 = 100 * q2,
        Location = step
      ) %>%
      mutate(Data = "Clinical trial observations")
  } else {
    dt_observed <-
      read_csv(file.path(dirX_inp(step), "data.csv")) %>%
      group_by(Drug_regimen, Day) %>%
      summarise(
        m  = 100 * binom.test(nr_recurrences, nr_enrolled)$estimate,
        q1 = 100 * binom.test(nr_recurrences, nr_enrolled)$conf.int[1],
        q2 = 100 * binom.test(nr_recurrences, nr_enrolled)$conf.int[2],
        Location = step
      ) %>%
      ungroup() %>%
      mutate(Data = "Clinical trial observations")
  }
  
  
  # prepare simulation results
  dt_simulated <-
    read_estimates(step, "MCMC", taskname) %>%
    handle_burnin(action = "remove") %>%
    add_known_efficacies(step, taskname) %>%
    transform_for_Kaplan_Meier() %>%
    add_input_variables(step, taskname) %>%
    rowwise() %>%
    group_split() %>%
    lapply(function (dt) {
      
      dt %>%
        # actual days from treatment onset
        bind_cols(Day = seq(from = 0, to = dt$maxday, by = 0.5)) %>%
        # remove BS prophylaxis period
        mutate(delta_t = Day - proph_BS) %>%
        # compute proportion of patients with BS by time t
        mutate(prop_I = if_else(delta_t < 1, 0, 100 * ll_bigF(delta_t, h_tP, gamma_L, rbite, rrela, 1 - eff * p_meta))) %>%
        # clean
        select(any_of(c("Location", "Drug_regimen", "Day", "prop_I")))
      
    }) %>%
    bind_rows() %>%
    group_by(pick(any_of(c("Location", "Drug_regimen", "Day")))) %>%
    summarise(
      m  = quantile(prop_I, probs = 0.5),
      q1 = quantile(prop_I, probs = 0.025),
      q2 = quantile(prop_I, probs = 0.975)
    ) %>%
    ungroup() %>%
    mutate(Data = "Model simulations")
  
  # combine
  dt <-
    bind_rows(dt_observed, dt_simulated) %>%
    mutate(Drug_regimen = factor(Drug_regimen, levels = label_drug_regimen$name))
  
  return(dt)
  
}


fig_Kaplan_Meier <- function (step, taskname, filter_dr = NULL, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": Kaplan-Meier ", taskname, "\n"))
  
  fig_name <- paste(step, "Kaplan_Meier", taskname, sep = "-")
  
  
  if (step %in% LETTERS) { # one plot for each step
    
    # Data for KM plots
    dt <- 
      prep_Kaplan_Meier_data(step, taskname) %>%
      mutate(drug_regimen = Drug_regimen)
    
    text_size <- 8
    ln_width <- 0.5
    
    if(step == "A") {
      y_lab <- expression(paste("Cumulative risk of symptomatic ", italic("P. vivax"), " recurrence (%)"))
    } else {
      y_lab <- expression(paste("Cumulative risk of ", italic("P. vivax"), " recurrence (%)"))
    }
    
  } else { # one plot combining all steps (poster)
    
    dt <-
      bind_rows(
        prep_Kaplan_Meier_data("A", taskname) %>% filter(grepl("PQ_high", Drug_regimen) & Location == "VN001"),
        prep_Kaplan_Meier_data("B", taskname) %>% filter(grepl("PQ_low",  Drug_regimen)),
        prep_Kaplan_Meier_data("C", taskname) %>% filter(grepl("TQ_",     Drug_regimen))
      ) %>%
      mutate(
        drug_regimen = 
          Drug_regimen %>%
          forcats::fct_rev() %>%
          {
            if (step == "poster") {
              relevel(., "PQ_highdose7") %>%
                relevel("PQ_highdose14")
            } else {
              relevel(., "PQ_lowdose")
            }
          } %>%
          forcats::fct_rev(),
        Location = if_else(Location == "VN001", "A", "BC")
      )
    
    text_size <- if (step == "poster") 21 else 25
    ln_width <- if (step == "talk_mini") 2 else 1
    pt_size <- if (step == "talk_mini") 7 else 4
    y_lab <- expression(paste(italic("P. vivax"), " recurrence (%)"))
    
  }
  
  
  # subset
  if (! is.null(filter_dr) ) {
    dt <- dt %>% filter(drug_regimen %in% filter_dr)
  }
  
  
  # plot Kaplan-Meier curves
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = Day, colour = Data)
    ) +
    # survival curves
    geom_step(mapping = aes(y = m),  linewidth = ln_width) +
    geom_line(mapping = aes(y = q1), linewidth = ln_width, linetype = "dashed") +
    geom_line(mapping = aes(y = q2), linewidth = ln_width, linetype = "dashed") +
    # scales
    scale_x_continuous(
      name = "Days from trial enrolment",
      minor_breaks = NULL
    ) +
    scale_y_continuous(name = y_lab) +
    theme_global(base_size = text_size) +
    theme(legend.position = if_else(step == "talk", "right", "bottom"))
  
  
  if (step == "A") {
    
    p <- p + facet_grid(rows = vars(Drug_regimen), cols = vars(Location), labeller = global_labeller)
    
    fig_width  <- 7.5
    fig_height <- 4
    
  } else if (step %in% c("B", "C")) {
    
    p <- 
      p +
      # error bars
      geom_errorbar(data = dt %>% filter(grepl("trial", Data)), mapping = aes(ymin = q1, ymax = q2), width = 5) +
      geom_point(data = dt %>% filter(grepl("trial", Data)), mapping = aes(y = m)) +
      # facets
      facet_grid(cols = vars(drug_regimen), labeller = global_labeller)
    
    fig_width  <- 7.5
    fig_height <- 2.5
    
  } else {
    
    p <- 
      p +
      # error bars
      geom_errorbar(data = dt %>% filter(grepl("trial", Data) & Location == "BC"), mapping = aes(ymin = q1, ymax = q2), linewidth = ln_width, width = 8) +
      geom_point(   data = dt %>% filter(grepl("trial", Data) & Location == "BC"), mapping = aes(y = m), size = pt_size) +
      # facets
      facet_wrap(
        vars(drug_regimen),
        labeller = global_labeller,
        scales = "free_x",
        dir = if_else(step == "poster", "h", "v"),
        as.table = if_else(step == "poster", FALSE, TRUE),
        ncol = 3
      ) +
      theme(legend.key.width = unit(35, "points"))
    
    
    fig_width  <- if_else(step == "poster", 10, 15)
    fig_height <- 7
    
  }
  
  
  # save to file
  my.ggsave(plot = p,
            path = dir_fig,
            filename = fig_name,
            device   = device,
            width    = fig_width,
            height   = fig_height)
  
}


fig_posterior_density <- function (step, taskname, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": posterior densities ", taskname, "\n"))
  
  fig_name <- paste(step, "posterior_densities", taskname, sep = "-")
  
  
  # read in
  if (step != "poster") {
    
    # one plot for each step
    dt <-
      read_estimates(step, "MCMC", taskname) %>%
      pivot_and_rescale() %>%
      compute_density()
    
    ln_width <- 0.75
    lab_title <- unique(dt$label_title)
    fig_width  <- 7.5
    fig_height <- ifelse(step == "A", 8, ifelse(step == "B", 2.5, 4))
    text_size <- 8
    
  } else {
    
    # one plot combining all steps (poster)
    dt <-
      c("B", "A", "C") %>%
      lapply(function (st) {
        
        read_estimates(st, "MCMC", taskname) %>%
          pivot_and_rescale() %>%
          compute_density() %>%
          # clean
          filter(grepl("eff__", name)) %>%
          select(x, y, m, y_median, position, label_facet2)
        
      }) %>%
      bind_rows() %>%
      rename(label_facet = label_facet2)
    
    ln_width <- 1.5
    lab_title <- NULL
    fig_width  <- 15
    fig_height <- 4
    text_size <- 20
    
  }
  
  
  # plot
  p <-
    ggplot() +
    ## lower tail
    geom_area(
      data = dt %>% filter(grepl("Lower", position)),
      mapping = aes(x = x, y = y, fill = position)
    ) +
    ## upper tail
    geom_area(
      data = dt %>% filter(grepl("Upper", position)),
      mapping = aes(x = x, y = y, fill = position)
    ) +
    ## middle 95%
    geom_area(
      data = dt %>% filter(grepl("Middle", position)),
      mapping = aes(x = x, y = y, fill = position)
    ) +
    ## median
    geom_segment(
      data = dt,
      mapping = aes(x = m, xend = m, y = y_median, yend = 0),
      linewidth = ln_width
    ) +
    facet_wrap(vars(label_facet), scales = "free", ncol = ifelse(step %in% c("C", "poster"), 3, 4)) +
    labs(title = lab_title) +
    scale_x_continuous(name = NULL, labels = number) +
    scale_y_continuous(name = "Posterior density", breaks = NULL) +
    scale_fill_discrete(
      name = "Probability",
      type = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0")
    ) +
    theme_global(base_size = text_size) +
    theme(legend.position = "bottom")
  
  
  # print to file
  my.ggsave(plot = p,
            path = dir_fig,
            filename = fig_name,
            device   = device,
            width    = fig_width,
            height   = fig_height)
  
}


fig_sensitivity_fits <- function (device = device_SI) {
  
  
  cat(paste0("# Fig A, B, C: sensitivity analysis\n"))
  
  fig_names <- paste("ABC-sensitivity", c("countries", "8AQ_action", "varying"), sep = "-")
  
  
  # read in
  dt <-
    tibble(step = c("A", "B", "C")) %>%
    # get tasknames for each step
    group_split(step) %>%
    lapply(function (row) {
      
      read_csv(file.path(dirX_inp(row$step), "varying.csv")) %>%
        mutate(taskname = factor(taskname, levels = taskname)) %>%
        bind_cols(row)
      
    }) %>%
    bind_rows() %>%
    # label names
    mutate(
      value = 
        case_when(
          grepl("CYP2D6", taskname) ~ paste("CYP2D6 prevalence =", CYP2D6_prev),
          grepl("dg",     taskname) ~ paste("Hypnozoite carriage =", d_gamma, "days"),
          grepl("CQ",     taskname) ~ paste("CQ prophylaxis =", proph_CQ, "days"),
          grepl("DP",     taskname) ~ paste("DP prophylaxis =", proph_DP, "days"),
          grepl("TQl",    taskname) ~ paste("Low-dose TQ prophylaxis =", proph_TQl, "days"),
          grepl("TQh",    taskname) ~ paste("High-dose TQ prophylaxis =", proph_TQh, "days"),
          .default                  = taskname
        )
    ) %>%
    select(step, taskname, value) %>%
    # compute summary for each task
    group_split(step, taskname, value) %>%
    lapply(function (row) {
      
      # read in
      read_estimates(row$step, "MCMC", row$taskname) %>%
        handle_burnin(action = "remove") %>%
        # select efficacy columns
        select(starts_with("eff")) %>%
        pivot_longer(
          cols = starts_with("eff"),
          names_prefix = "eff__",
          names_to = "drug_regimen",
          values_to = "eff"
        ) %>%
        mutate(eff = 100 * eff) %>%
        # pivot_and_rescale() %>%
        bind_cols(row)
      
    }) %>%
    bind_rows() %>%
    # compute summary
    group_by(step, taskname, value, drug_regimen) %>%
    summarise(
      m  = median(eff),
      q1 = quantile(eff, probs = 0.025),
      q2 = quantile(eff, probs = 0.975)
    ) %>%
    ungroup() %>%
    # clean
    mutate(
      step = factor(step),
      drug_regimen = factor(drug_regimen, label_drug_regimen$name)
    )
  
  
  # multiple plots
  for (fig_name in fig_names) {
    
    
    # filter and rename
    if (grepl("countr", fig_name)) {
      
      dt_plot <-
        dt %>%
        filter(
          step == "A" &
            taskname %in% c("AON", "Afghanistan", "Ethiopia", "Indonesia", "Vietnam")
        ) %>%
        mutate(
          value = if_else(value == "AON", "All countries *", value)
        )
      
    } else if (grepl("action", fig_name)) {
      
      dt_plot <-
        dt %>%
        filter(taskname %in% c("AON", "leaky")) %>%
        mutate(value = if_else(value == "AON", "AON *", value))
      
    } else {
      
      dt_plot <-
        dt %>%
        filter(grepl("AON|minus|plus", taskname)) %>%
        mutate(value = if_else(value == "AON", "Default values *", value))
      
    }
    
    
    # plot
    p <- 
      ggplot(
        data = dt_plot,
        mapping = aes(x = taskname, y = m, ymin = q1, ymax = q2)
      ) +
      geom_pointrange() +
      scale_y_continuous(
        name = "Hypnozoiticidal efficacy (%)",
        limits = c(0, 100)
      ) +
      scale_x_discrete(
        name = NULL,
        breaks = dt_plot$taskname,
        labels = dt_plot$value
      ) +
      facet_wrap(
        vars(drug_regimen),
        labeller = global_labeller,
        scales = "free_x",
        nrow = 1
      ) +
      theme_global(base_size = 8)
    
    
    if (grepl("varying", fig_name)) {
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    
    # print to file
    my.ggsave(plot = p,
              path = dir_fig,
              filename = fig_name,
              device   = device,
              width    = 7.5,
              height   = 3) 
    
  }
  
}


prop_aged_less <- function (age_in_days, age_mean, age_max) {
  
  (1 - exp(- age_in_days / age_mean)) / (1 - exp(- age_max / age_mean))
  
}


prop_G6PD_activity_less <- function (sex, threshold, G6PD_prev, mu_G6PD_nor, sig_G6PD_nor, mu_G6PD_het, sig_G6PD_het, mu_G6PD_def, sig_G6PD_def) {
  
  if (sex == "male") {
    
    # hemizygous deficient male
    G6PD_prev * pnorm(threshold, mean = mu_G6PD_def, sd = sig_G6PD_def) +
      # hemizygous normal male
      (1 - G6PD_prev) * pnorm(threshold, mean = mu_G6PD_nor, sd = sig_G6PD_nor)
    
  } else if (sex == "female") {
    
    # homozygous deficient female
    G6PD_prev^2 * pnorm(threshold, mean = mu_G6PD_def, sd = sig_G6PD_def) +
      # heterozygous deficient female
      2 * G6PD_prev * (1 - G6PD_prev) * pnorm(threshold, mean = mu_G6PD_het, sd = sig_G6PD_het) +
      # homozygous normal female
      (1 - G6PD_prev)^2 * pnorm(threshold, mean = mu_G6PD_nor, sd = sig_G6PD_nor)
    
  } else {
    NA
  }
  
}


fig_survival_curves <- function (step, taskname, filter_dr = NULL, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": survival curves\n"))
  
  fig_name <- paste(step, "survival_curves", sep = "-")
  
  
  if (step == "ABC") {
    y_lab <- expression(paste("Cumulative risk of ", italic("P. vivax"), " recurrence (%)"))
    text_size <- 8
    ln_width <- 0.6
    legend_rows <- 2
    fig_width  <- 7.5
    fig_height <- 4
  } else {
    y_lab <- expression(paste(italic("P. vivax"), " recurrence (%)"))
    text_size <- 25
    ln_width <- 1.5
    legend_rows <- if (step == "talk") 2 else 1
    fig_width  <- 15
    fig_height <- 8.5
  }
  
  # fix survival function
  ll_bigF <- if (taskname == "leaky") leaky_bigF else AON_bigF
  
  # get estimates
  file_path <- file.path(dir_tab, "ABC-estimates-eff-AON.csv")
  if (! file.exists(file_path))
    tab_estimates_efficacy_effectiveness(taskname = "AON")
  
  dt_eff <- 
    read_csv(file_path) %>%
    filter(group != "Description") %>%
    select(contains("_"))
  
  el <- read_csv(file.path(dir_tab, "proportion_eligible_for_8AQ.csv"))
  
  etc <-
    dt_eff %>%
    select(row_, etc_m) %>%
    pivot_wider(names_from = row_, values_from = etc_m, values_fn = ~ .x / 100)
  
  
  dt <-
    dt_default %>%
    # expand tibble
    expand_grid(
      # vary transmission intensity
      yEIR2 = c(0.1, 1, 10),
      # actual days from treatment onset
      Day = seq(from = 0, to = 180, by = 0.5)
    ) %>%
    # compute helper values
    mutate(
      # daily rates
      gamma  = 1 / d_gamma,
      lambda = yEIR2 / 365,
      ff     = 1 / d_ff,
      # hazards
      h_tP_CQ  = exp( - gamma * CQ__proph),
      h_tP_TQl = exp( - gamma * TQ_lowdose__proph),
      h_tP_TQh = exp( - gamma * TQ_highdose__proph),
      # remove prophylaxis period
      Day_CQ  = Day - CQ__proph,
      Day_TQl = Day - TQ_lowdose__proph,
      Day_TQh = Day - TQ_highdose__proph
    ) %>%
    mutate(
      
      # trial conditions: compute proportion of patients with BS infection by time t
      No_8AQs__tc           = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$No_8AQs) ),
      AQ_ideal__tc          = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$AQ_ideal) ),
      
      TQ_lowdose__tc        = if_else( Day_TQl < 1, 0, ll_bigF(Day_TQl, h_tP_TQl, gamma, lambda, ff, 1 - etc$TQ_lowdose) ),
      TQ_highdose__tc       = if_else( Day_TQh < 1, 0, ll_bigF(Day_TQh, h_tP_TQh, gamma, lambda, ff, 1 - etc$TQ_highdose) ),
      
      PQ_lowdose14__tc      = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$PQ_lowdose14) ),
      PQ_lowdose7__tc       = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$PQ_lowdose7) ),
      PQ_highdose14__tc     = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$PQ_highdose14) ),
      PQrestr_highdose7__tc = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - etc$PQrestr_highdose7) ),
      
      # temporary: PQ accounting for PQ adherence
      PQ_lowdose14__oc_temp     = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - (etc$PQ_lowdose14      * PQ_14__adhere)) ),
      PQ_lowdose7__oc_temp      = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - (etc$PQ_lowdose7       * PQ_7__adhere )) ),
      PQ_highdose14__oc_temp    = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - (etc$PQ_highdose14     * PQ_14__adhere)) ),
      PQrestr_highdose__oc_temp = if_else( Day_CQ  < 1, 0, ll_bigF(Day_CQ,  h_tP_CQ,  gamma, lambda, ff, 1 - (etc$PQrestr_highdose7 * PQ_7__adhere )) ),
      
      # operational conditions: compute proportion of patients with BS infection by time t
      No_8AQs__oc           = No_8AQs__tc,
      AQ_ideal__oc          = AQ_ideal__tc,
      
      PQ_lowdose14__oc      = (1 - el$elig_PQ) * No_8AQs__oc + el$elig_PQ * PQ_lowdose14__oc_temp,
      PQ_lowdose7__oc       = (1 - el$elig_PQ) * No_8AQs__oc + el$elig_PQ * PQ_lowdose7__oc_temp,
      PQ_highdose14__oc     = (1 - el$elig_PQ) * No_8AQs__oc + el$elig_PQ * PQ_highdose14__oc_temp,
      
      PQrestr_highdose7__oc = (1 - el$elig_PQ) * No_8AQs__oc + (el$elig_PQ - el$elig_PQrestr) * PQ_lowdose7__oc_temp + el$elig_PQrestr * PQrestr_highdose__oc_temp,
      TQ_lowdose__oc        = (1 - el$elig_PQ) * No_8AQs__oc + (el$elig_PQ - el$elig_TQ     ) * PQ_lowdose7__oc_temp + el$elig_TQ      * TQ_lowdose__tc,
      TQ_highdose__oc       = (1 - el$elig_PQ) * No_8AQs__oc + (el$elig_PQ - el$elig_TQ     ) * PQ_lowdose7__oc_temp + el$elig_TQ      * TQ_highdose__tc
      
      
    ) %>%
    # clean
    select(Day, yEIR2, ends_with(c("__tc", "__oc"))) %>%
    pivot_longer(
      cols = matches("__"),
      names_to = c("drug_regimen", "conditions"),
      names_sep = "__"
    ) %>%
    mutate(
      conditions   = factor(conditions,   levels = names(label_conditions)),
      Conditions   = factor(conditions,   levels = names(label_conditions2)),
      drug_regimen = factor(drug_regimen, levels = rev(label_drug_regimen$name)),
      yeir2 = yEIR2,
      value = 100 * value
    )
  
  
  # prophylaxis table for plotting vertical dashed lines
  dt_proph <-
    tibble(
      drug_regimen = c("No_8AQs", "TQ_lowdose", "TQ_highdose"),
      proph        = c(
        dt_default$CQ__proph,
        dt_default$TQ_lowdose__proph,
        dt_default$TQ_highdose__proph
      )
    )
  
  
  # subset
  if (! is.null(filter_dr) ) {
    dt       <- dt       %>% filter(drug_regimen %in% filter_dr)
    dt_proph <- dt_proph %>% filter(drug_regimen %in% filter_dr)
  } else {
    dt <- dt %>% filter( !(drug_regimen %in% c("AQ_ideal")))
  }
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = Day, y = value, colour = drug_regimen)
    ) +
    # duration of post-treatment prophylaxis
    geom_vline(
      data = dt_proph,
      mapping = aes(xintercept = proph, colour = drug_regimen),
      linetype = "dashed",
      linewidth = ln_width,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    # cumulative incidence
    geom_line(linewidth = ln_width) +

    labs(x = "Days from start of treatment", y = y_lab) +
    my.scale_drug_regimen(aesthetics = "colour", na.value = "black") +
    guides(colour = guide_legend(reverse = TRUE, nrow = legend_rows)) +
    theme_global(base_size = text_size) +
    theme(legend.position = "bottom")
  
  
  if (step == "ABC") {
    p <- p + 
      facet_grid(rows = vars(conditions), cols = vars(yEIR2), labeller = global_labeller)
  } else {
    p <- p + 
      facet_grid(rows = vars(Conditions), cols = vars(yeir2), labeller = global_labeller) +
      theme(legend.key.size = unit(32, "points"))
    
  }
  
  
  # print to file
  my.ggsave(plot = p,
            path = dir_fig,
            filename = fig_name,
            device   = device,
            width    = fig_width,
            height   = fig_height)
  
}
