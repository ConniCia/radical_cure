read_in <- function (step, file_name) {
  
  # read in
  dt <-
    read_csv(
      file = file.path(dirX_res(step), file_name),
      col_select = -starts_with("path_")
    )
  
  
  # transform
  if (file_name == "effsize.csv") {
    
    dt <-
      dt %>%
      # filter
      filter(-0.5 < year) %>%
      # cumulative sum of cases
      group_by(across(-any_of(c("year", "sum_new_T", "mean_prev_PCR")))) %>%
      arrange(year) %>%
      mutate(cumsum_new_T = cumsum(sum_new_T)) %>%
      ungroup() %>%
      select(-sum_new_T)
    
    dt1 <- dt %>% filter(drug_regimen != "No_8AQs")
    dt2 <-
      dt %>%
      filter(drug_regimen == "No_8AQs") %>%
      select(-drug_regimen, -pars_varying, -PQ_7__adhere, -PQ_14__adhere, -CYP2D6_prev, -G6PD_prev) %>%
      rename(
        baseline_mean_prev_PCR = mean_prev_PCR,
        baseline_cumsum_new_T  = cumsum_new_T
      )
    
    dt <-
      # compute reductions w.r.t CQ-only drug regimen
      inner_join(dt1, dt2, by = dt2 %>% select(-starts_with("baseline_")) %>% names()) %>%
      mutate(
        effsize       =  100 * (baseline_mean_prev_PCR - mean_prev_PCR) / baseline_mean_prev_PCR,
        cases_averted = 1000 * (baseline_cumsum_new_T  - cumsum_new_T)  / N_pop,
      ) %>%
      select(-contains(c("sample_nr", "mean_prev_PCR", "cumsum_new_T"))) %>%
      # pivot longer
      pivot_longer(cols = c(effsize, cases_averted)) %>%
      # compute summary stats: median and 95% CrI
      group_by(across(-"value")) %>%
      summarise(
        m  = median  (value,                na.rm = TRUE),
        q1 = quantile(value, probs = 0.025, na.rm = TRUE),
        q2 = quantile(value, probs = 0.975, na.rm = TRUE)
      ) %>%
      ungroup()
    
  } else {
    
    # decompress results
    dt <-
      dt %>%
      pivot_longer(
        cols = matches("m__|q1__|q2__"),
        names_to = c("stat", "variable"),
        names_sep = "__",
        values_to = "value"
      ) %>%
      pivot_wider(
        names_from  = "stat",
        values_from = "value"
      ) %>%
      # only save non-zero results
      filter(q2 > 0.5)
    
  }
  
  
  # clean
  dt <-
    dt %>%
    mutate(
      # define drug regimen factor levels
      drug_regimen = factor(drug_regimen, label_drug_regimen$name),
      # add seasonality labels
      seasonality = cut(
        kappa_seas,
        breaks = c(0, 0.9, 1.1, 100),
        labels = c("None", "Moderate", "High")
      ),
      # replace PQ adherence labels
      PQ__adhere = cut(
        PQ_7__adhere,
        breaks = c(0, 0.5, 0.7, 1),
        labels = c("Low", "Moderate", "High")
      ),
      # PQ_14__adhere = NULL
    )
  
}


fig_trends <- function (dt_timeline, year_endpoint, filter_seas, device = device_SI) {
  
  
  cat(paste0("# Fig D: trends, ", filter_seas, " seasonality\n"))
  
  # figure name
  fig_name <- paste("D-trends", filter_seas, sep = "-")
  
  
  # transform
  dt <-
    dt_timeline %>%
    # filter
    filter(
      -3 < year &
        d_ff == dt_default$d_ff &
        PQ__adhere == "Moderate" &
        seasonality == filter_seas &
        grepl("default|yEIR|seas", pars_varying)
    ) %>%
    # transform
    mutate(
      m =
        if_else(
          grepl("new_T", variable),
          365 * 1000 * m / N_pop,
          100 * m / N_pop
        )
    ) %>%
    select(year, drug_regimen, yEIR, variable, m) %>%
    mutate(drug_regimen = forcats::fct_rev(drug_regimen))
  
  
  # plot
  p <-
    ggplot(
      data = dt %>% arrange(drug_regimen),
      mapping = aes(x = year, y = m, colour = drug_regimen)
    ) +
    # trend
    geom_line() +
    # start of intervention
    geom_vline(xintercept = dt_default$intervention_time) +
    geom_vline(xintercept = dt_default$intervention_time + year_endpoint, linetype = "dashed") +
    # facet
    facet_grid(
      rows = vars(variable),
      cols = vars(yEIR),
      labeller = global_labeller,
      scales = "free",
      switch = "y"
    ) +
    # layout
    scale_x_continuous(
      name = "Years from introduction of radical cure in case management",
      breaks = seq(-2, 10, 2)
    ) +
    scale_y_continuous(name = NULL, limits = c(0, NA)) + # always include 0 in y-axis
    my.scale_drug_regimen(aesthetics = "colour") +
    guides(colour = guide_legend(reverse = TRUE)) +
    theme_global(base_size = 8) +
    theme(
      legend.position = "bottom",
      strip.placement = "outside",
      strip.background.y = element_blank()
    )
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = 7.5,
            height   = 3.5)
  
}


fig_effsize <- function (step, dt_effsize, year_endpoint, filter_dr = NULL, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": effect size\n"))
  
  tab_name <- paste0(step, "-cite-effsize")
  fig_name <- paste0(step, "-effsize")
  
  
  # data
  dt <-
    dt_effsize %>%
    # filter
    filter(
      year == dt_default$intervention_time + year_endpoint &
        d_ff == dt_default$d_ff &
        seasonality == "Moderate" &
        !(grepl("TQ", drug_regimen) & PQ__adhere != "Moderate") &
        grepl("default|yEIR|PQ__adhere", pars_varying) &
        name == "effsize"
    ) %>%
    # transform
    mutate(
      yeir = yEIR,
      PQ__adhere = factor(PQ__adhere, c("High", "Moderate", "Low"))
    ) %>%
    select(m, q1, q2, drug_regimen, yEIR, yeir, PQ__adhere)
  
  
  # if (step != "D") {
  #   dt <- dt %>% filter(PQ__adhere == "Moderate")
  # }
  
  
  # save table for citation
  if (step == "D") {
    dt %>%
      select(drug_regimen, yEIR, PQ__adhere, m, q1, q2) %>%
      mutate(
        CrI = paste0(round(m, 1), "% (", round(q1, 1), "-", round(q2, 1), ")"),
        m   = round(m),
        q1  = NULL,
        q2  = NULL
      ) %>%
      # print to file
      write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  }
  
  
  # subset
  if (! is.null(filter_dr) ) {
    dt <- dt %>% filter(drug_regimen %in% filter_dr)
  }
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      aes(x = m, xmin = q1, xmax = q2,
          y = drug_regimen, fill = drug_regimen,
          group = PQ__adhere,
          pattern = PQ__adhere, pattern_angle = PQ__adhere)
    ) +
    # column plot
    geom_col_pattern(
      position = position_dodge(preserve = "total"),
      colour = "white", 
      pattern_colour = "grey80",
      pattern_density = if (step == "D") 0.1   else 0.2,
      pattern_spacing = if (step == "D") 0.025 else 0.03,
      pattern_key_scale_factor = 0.6
    ) + 
    geom_errorbar(position = position_dodge(width = 0.9), colour = "black", width = 0.5, linewidth = if_else(step == "D", 0.5, 1.5)) +
    # layout
    coord_cartesian(xlim = c(0, 100)) +
    scale_x_continuous(
      name = paste("Reduction in PCR-prevalence after", year_endpoint, "years (%)"),
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    my.scale_drug_regimen(aesthetics = "fill") +
    scale_pattern_manual(
      name   = "PQ adherence",
      values = c(Low = "stripe", Moderate = "none", High = "stripe")
    ) +
    scale_pattern_angle_manual(
      name   = "PQ adherence",
      values = c(Low = -45, Moderate = 0, High = 45)
    ) +
    guides(
      pattern = guide_legend(
        override.aes = list(fill = "white", colour = "black", pattern_colour = "grey70")
      ),
      fill = guide_legend(override.aes = list(pattern = "none"))
    ) +
    theme_global(base_size = if_else(step == "D", 8, 23)) +
    theme(legend.position = "bottom")
  
  
  if (step == "D") {
    
    p <-
      p +
      # facet
      facet_wrap(vars(yEIR), labeller = global_labeller) +
      # layout
      scale_y_discrete(
        name = NULL,
        limits = rev, # reverse the y axis but not the order in the legend
        breaks = NULL
      ) +
      theme(
        legend.key.size = unit(10, "points"),
        legend.box = "vertical",
        # legend.spacing = unit(40, "points"),
        # legend.key = element_rect(size = 100),
        # legend.box.just = "left"
      )
    
    fig_width <- 7.5
    fig_height <- 3.5
    
  } else {
    
    p <-
      p +
      # facet
      facet_wrap(vars(yeir), labeller = global_labeller) +
      # layout
      scale_y_discrete(
        name = NULL,
        limits = rev, # reverse the y axis but not the order in the legend
        labels = tibble::deframe(label_drug_regimen %>% select(name, label))
      ) +
      guides(fill = "none") +
      theme(
        legend.key.size = unit(30, "points"),
        axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank()
      )
    
    fig_width <- 16
    fig_height <- if (step == "poster") 6 else 7
    
  }
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = fig_width,
            height   = fig_height)
  
}


fig_cases_averted <- function (dt_effsize, year_endpoint, device = device_SI) {
  
  
  cat(paste0("# Fig D: cases averted\n"))
  
  tab_name <- "D-cite-cases_averted"
  fig_name <- "D-cases_averted"
  
  
  # data
  dt <-
    dt_effsize %>%
    # filter
    filter(
      pars_varying %in% c("default", "yEIR") &
        name == "cases_averted"
    ) %>%
    # transform
    select(year, m, q1, q2, drug_regimen, yEIR)
  
  
  # save table for citation
  dt %>%
    filter(year == year_endpoint) %>%
    select(drug_regimen, yEIR, m) %>%
    mutate(m = round(m)) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = year, y = m, colour = drug_regimen)
    ) +
    # cases averted
    geom_line(linewidth = 0.6) +
    # vertical line
    geom_vline(xintercept = dt_default$intervention_time + year_endpoint, linetype = "dashed") +
    # facet
    facet_wrap(vars(yEIR), labeller = global_labeller) +
    # layout
    scale_x_continuous(
      name = "Years from change in case management",
      breaks = seq(0, 10, 2)
    ) +
    scale_y_continuous(
      name = "Cumulative averted clinical cases\n(per 1 000 population)",
      labels = number
    ) +
    my.scale_drug_regimen(aesthetics = "colour") +
    theme_global(base_size = 8) +
    theme(legend.position = "bottom")
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = 7.5,
            height   = 2.5)
  
}


fig_treatment_courses <- function (dt_stats, device = device_SI) {
  
  
  cat(paste0("# Fig D: treatment courses\n"))
  
  fig_name <- "D-treatment_courses"
  
  
  # data
  dt <-
    dt_stats %>%
    # filter
    filter(
      -0.5 < year &
        variable %in% label_treatment_courses$name &
        pars_varying %in% c("default", "yEIR")
    ) %>%
    # transform
    mutate(
      # rescale
      m = 1000 * m / N_pop,
      # special case AQ_ideal
      variable = if_else(drug_regimen == "AQ_ideal" & variable == "sum_new_PQ", "sum_new_AQ", variable)
    ) %>%
    # special case PQrestr_highdose7
    group_by(drug_regimen, year, yEIR) %>%
    mutate(
      variable = if_else(drug_regimen == "PQrestr_highdose7" & variable == "sum_new_TQ", "sum_new_PQ", variable)
    ) %>%
    ungroup() %>%
    # clean
    mutate(variable = factor(variable, label_treatment_courses$name)) %>%
    rename(Drug_regimen = drug_regimen) %>%
    select(year, Drug_regimen, yEIR, variable, m)
  
  
  # plot
  p <-
    ggplot() +
    # background
    geom_rect(
      data = tibble(a = seq(0, 10, 2) - 0.5, b = seq(0, 10, 2) + 0.5),
      aes(xmin = a, xmax = b, ymin = -Inf, ymax = Inf),
      fill = "grey60", alpha = 0.4
    ) +
    # effect size
    geom_bar(
      data = dt,
      mapping = aes(x = year, weight = m, fill = variable),
      position = position_dodge(preserve = "single")
    ) +
    # facet
    facet_grid(
      rows = vars(Drug_regimen),
      cols = vars(yEIR),
      labeller = global_labeller
    ) +
    # layout
    geom_hline(yintercept = Inf, color = "white", linewidth = 1) +
    scale_x_continuous(name = "Years from change in case management", breaks = 0:10) +
    scale_y_continuous(
      name = "Number of tests or treatment courses (per year per 1 000 population)",
      labels = number,
      limits = c(0, 580)
    ) +
    my.scale_treatment_courses(aesthetics = "fill") +
    theme_global(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(10, "points"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = 7.5,
            height   = 8.5)
  
}


fig_relapses_reinfections <- function (dt_stats, device = device_SI) {
  
  
  cat(paste0("# Fig D: relapses and reinfections\n"))
  
  fig_name <- "D-relapses_reinfections"
  
  
  # data
  dt <-
    dt_stats %>%
    # filter
    filter(
      -0.5 < year &
        variable %in% c("sum_new_CQ_bite", "sum_new_CQ_relapse") &
        pars_varying %in% c("default", "yEIR")
    ) %>%
    # transform
    mutate(
      # rescale
      m = 1000 * m / N_pop
    ) %>%
    # clean
    mutate(variable = factor(variable, label_relapses_reinfections$name)) %>%
    rename(Drug_regimen = drug_regimen) %>%
    select(year, Drug_regimen, yEIR, variable, m)
  
  
  # plot
  p <-
    ggplot() +
    # background
    geom_rect(
      data = tibble(a = seq(0, 10, 2) - 0.5, b = seq(0, 10, 2) + 0.5),
      aes(xmin = a, xmax = b, ymin = -Inf, ymax = Inf),
      fill = "grey60", alpha = 0.4
    ) +
    # effect size
    geom_col(
      data = dt,
      mapping = aes(x = year, y = m, fill = variable),
      position = "dodge"
    ) +
    # facet
    facet_grid(
      rows = vars(Drug_regimen),
      cols = vars(yEIR),
      labeller = global_labeller
    ) +
    # layout
    geom_hline(yintercept = Inf, color = "white", linewidth = 1) +
    scale_x_continuous(name = "Years from change in case management", breaks = 0:10) +
    scale_y_continuous(
      name = "Number of clinical cases (per year per 1000 population)",
      labels = number,
      limits = c(0, 450)
    ) +
    my.scale_relapses_reinfections(aesthetics = "fill") +
    theme_global(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(10, "points"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = 7.5,
            height   = 8.5)
  
}


fig_seasonality <- function (dt_timeline, device = device_SI) {
  
  
  cat(paste0("# Fig D: seasonality patterns\n"))
  
  fig_name <- "D-seasonality"
  
  
  # read in and transform
  dt <-
    dt_timeline %>%
    select(dry_seas, kappa_seas, seasonality) %>%
    unique() %>%
    group_split(dry_seas, kappa_seas, seasonality) %>%
    lapply(function(x) {
      
      t <- seq(0, 3, 0.01) # time in years
      y <- x$dry_seas + (1 - x$dry_seas) * g(x$kappa_seas, t)
      
      bind_cols(x, t = t, y = y)
      
    }) %>%
    bind_rows() %>%
    rename(Seasonality = seasonality)
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = t, y = y, colour = Seasonality)
    ) +
    geom_line() +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mosquito density", breaks = NULL) +
    theme_global(base_size = 8)
  
  
  # save to file
  my.ggsave(plot     = p,
            path     = dir_fig,
            filename = fig_name,
            device   = device,
            width    = 7.5,
            height   = 3)
  
}


fig_sensitivity_IBM <- function (dt_effsize, dt_stats, year_endpoint, device = device_SI) {
  
  
  cat(paste0("# Fig D: sensitivity analyses\n"))
  
  fig_names <- paste("D-sensitivity", c("transmission_settings", "intervention_characteristics", "immunity_mechanisms"), sep = "-")
  
  
  # filter
  dt1 <-
    dt_effsize %>%
    filter(year == dt_default$intervention_time + year_endpoint)
  
  dt2 <-
    dt_stats %>%
    filter(
      year == dt_default$intervention_time &
        variable == "mean_prev_PCR" &
        (drug_regimen == "No_8AQs" | pars_varying %in% c("PQ__adhere", "CYP2D6_prev", "G6PD_prev"))
    ) %>%
    mutate(
      drug_regimen = NA,
      across(c(m, q1, q2), ~ 100 * .x / N_pop)
    ) %>%
    unique() %>%
    rename(name = variable)
  
  dt <- bind_rows(dt1, dt2) %>% filter(!grepl("-", pars_varying))
  
  
  # separate default and non-default parameter values
  dt_std <- dt %>% filter(pars_varying == "default") %>% mutate(type = "default")
  dt     <- dt %>% filter(pars_varying != "default") %>% mutate(type = "varying")
  
  
  # add default parameter values to each group of varying parameters
  dt <-
    dt %>%
    group_split(pars_varying) %>%
    lapply(function (x) {
      
      nv <- unique(x$pars_varying)
      
      x %>%
        # bind rows of jobs with default values
        bind_rows(dt_std %>% mutate(pars_varying = nv)) %>%
        # prep column of parameter whose value changes
        mutate(value = as.character(.data[[nv]])) %>%
        # what columns to keep  
        select(pars_varying, type, value, m, q1, q2, drug_regimen, name)
      
    }) %>%
    bind_rows() %>%
    # order values and denote default values
    mutate(
      value_num = case_when(
        value == "Low"      ~ -4,
        value == "None"     ~ -3,
        value == "Moderate" ~ -2,
        value == "High"     ~ -1,
        TRUE                ~ as.numeric(value)
      ),
      value_chr = if_else(type == "default", paste(value, "*"), value)
    ) %>%
    arrange(value_num) %>%
    # split into three figures
    mutate(
      fig = case_when(
        pars_varying %in% c("CM_cover", "PQ__adhere") ~ fig_names[2],
        grepl("_rel", pars_varying)                   ~ fig_names[3],
        TRUE                                          ~ fig_names[1]
      ),
      name         = factor(name, names(label_name)),
      pars_varying = factor(pars_varying, names(label_pars_varying))
    ) %>%
    # only show PQ drug regimens when varying PQ adherence
    filter(!(grepl("AQ|TQ", drug_regimen) & pars_varying == "PQ__adhere")) %>%
    # shaded background
    group_by(fig, pars_varying, value_num) %>%
    mutate(x = cur_group_id()) %>%
    ungroup()
  
  
  # add shaded background
  dt <-
    bind_rows(
      dt,
      dt %>%
        select(fig, pars_varying, name, value_num, value_chr, x) %>%
        unique() %>%
        arrange(x) %>%
        mutate(
          xa = if_else(x %% 2 == 1, x - 0.5, NA_real_),
          xb = if_else(x %% 2 == 1, x + 0.5, NA_real_)
        )
    )
  
  
  # plot and save to file
  for (fig_name in fig_names) {
    
    
    dt_plot <- dt %>% filter(fig == fig_name)
    
    # plot
    p <-
      ggplot(
        data = dt_plot,
        mapping = aes(x = x, y = m, xmin = xa, xmax = xb, ymin = q1, ymax = q2,
                      group = drug_regimen, colour = drug_regimen)
      ) +
      # shaded background
      geom_rect(
        ymin = -Inf, ymax = Inf,
        colour = NA, fill = "grey60", alpha = 0.4
      ) +
      # points
      geom_pointrange(data = dt_plot %>% filter(!is.na(drug_regimen)),
                      position = position_dodge(width = 1),
                      size = 0.15, linewidth = 0.5) +
      geom_pointrange(data = dt_plot %>% filter(is.na(drug_regimen)),
                      colour = "black",
                      position = position_dodge(width = 1),
                      size = 0.15, linewidth = 0.5) +
      # facets
      facet_grid(
        rows = vars(name),
        cols = vars(pars_varying),
        labeller = global_labeller,
        scales = "free",
        switch = "y"
      ) +
      scale_x_continuous(name = NULL, breaks = dt$x, labels = dt$value_chr) +
      scale_y_continuous(name = NULL, limits = c(0, NA), labels = number) + # always include 0 in y-axis
      geom_hline(yintercept = Inf, color = "white", linewidth = 1) +
      my.scale_drug_regimen(aesthetics = "colour") +
      theme_global(base_size = 8) +
      theme(
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside",
        strip.background.y = element_blank(),
        # panel.spacing.x = unit(4, "points"),
        # plot.margin   = margin(1,1,0,0)
      )
    
    
    # save to file
    my.ggsave(plot     = p,
              path     = dir_fig,
              filename = fig_name,
              device   = device,
              width    = 7.5,
              height   = 5)
    
  }  
  
}
