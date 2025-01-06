read_estimates <- function (step, filetype, taskname) {
  
  read_csv(file.path(dirX_res(step), paste0(filetype, "-", taskname, ".csv"))) %>%
    # clean
    mutate(
      i     = if ("i"     %in% colnames(.)) as.integer(i)    else row_number(),
      chain = if ("chain" %in% colnames(.)) as.factor(chain) else NULL
    )
  
}


compute_proportion_of_relapses <- function (dt) {
  
  for (lab in label_location$name) {
    
    if (any(grepl(lab, names(dt)))) {

      dt[[paste0("prela__", lab)]] <- dt[[paste0("rrela__", lab)]] / (dt[[paste0("rrela__", lab)]] + dt[[paste0("rbite__", lab)]])
      
    }
    
  }
  
  return(dt)
  
}


pivot_and_rescale <- function (dt) {
  
  dt %>%
    # pivot
    pivot_longer(cols = starts_with(c("rbite", "rrela", "prela", "eff", "ll"))) %>%
    # rescale parameter values
    mutate(
      value = 
        case_when(
          grepl("^eff",   name) ~ value * 100,
          grepl("^rbite", name) ~ value * 365,
          grepl("^rrela", name) ~ value * 365,
          grepl("^prela", name) ~ value * 100,
          TRUE                  ~ value
        )
    )
  
}


compute_summary_stats <- function (dt) {
  
  # remove burnin
  dt <- dt %>% handle_burnin(action = "remove")
  
  
  # compute multivariate ESS
  multi_ESS <-
    dt %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    select(-c(ll, i, chain)) %>%
    mcmcse::multiESS() %>%
    number(accuracy = 1)
  
  min_ESS <- mcmcse::minESS(p = length(unique(dt$name)) - 1) %>% number()
  
  
  # compute other summary stats
  dt %>%
    # merge chains
    select(-i, -chain) %>%
    # compute summary
    group_by(name) %>%
    summarise(
      m   = median(value),
      q1  = quantile(value, probs = 0.025),
      q2  = quantile(value, probs = 0.975),
      ESS = value %>% mcmcse::ess() %>% round()
    ) %>%
    ungroup() %>%
    # clean
    separate_wider_delim(
      name,
      delim = "__",
      names = c("temp1", "temp2"),
      too_few = "align_start",
      cols_remove = FALSE
    ) %>%
    mutate(
      # names
      line1 = label_fitted[temp1],
      line2 = c(
        tibble::deframe(label_location     %>% select(name, label2)),
        tibble::deframe(label_drug_regimen %>% select(name, label ))
        )[temp2],
      line2 = if_else(is.na(line2), "", line2),
      # CrI
      nd = if_else(q1 > 9, 1, 2),
      line3 = paste0(round(m, nd), " [", round(q1, nd), ", ", round(q2, nd), "]"),
      # ESS
      line4 = if_else(name == "ll", "", paste0("ESS = ", number(ESS, accuracy = 1))),
      # factors
      name2 = factor(name,  levels = order_fitted_pars),
      name3 = factor(temp2, levels = c(label_location$name, label_drug_regimen$name))
    ) %>%
    arrange(name2) %>%
    mutate(
      label_facet  = factor(name2, labels = paste(line1, line2, line3, line4, sep = "\n")),
      label_facet2 = factor(name2, labels = paste(line2, line3, sep = "\n")),
      label_title = paste0(
        "Multivariate ESS = ", multi_ESS, " (minimal multivariate ESS = ", min_ESS, ")"
      )
    ) %>%
    select(-nd, -starts_with("temp"))
  
}


thin_chains <- function (dt) {
  
  # at what interval should I sample if I wanted 500 points?
  delta_thin <- (max(dt$i) - min(dt$i)) / 500
  
  # round to the nearest 10s
  delta_thin <- as.integer(delta_thin / 10) * 10L
  
  # compute indices to be sampled
  i_thinned <- seq(from = min(dt$i), to = max(dt$i), by = delta_thin)
  
  
  dt %>% filter(i %in% i_thinned)
  
}


handle_burnin <- function (dt, action) {
  
  if (action == "remove") {
    dt %>% filter(i > burnin * max(i))
  } else if (action == "extract") {
    dt %>% filter(i < burnin * max(i))
  }
  
}


read_data_A <- function (step) {
  
  dt_observed <-
    read_csv(file.path(dirX_inp(step), "data.csv")) %>%
    mutate(across(where(is.numeric), as.integer)) %>%
    # compute nr_enrolled in each arm and site
    group_by(Location, Drug_regimen) %>%
    mutate(nr_enrolled = n()) %>%
    ungroup() %>%
    # count nr of dropouts and disease events at each time step
    group_by(Location, Drug_regimen, nr_enrolled, Day) %>%
    summarise(
      nr_disease = sum(fail),
      nr_dropout = n() - sum(fail)
    ) %>%
    ungroup() %>%
    # compute cumulative disease and dropout events
    group_by(Location, Drug_regimen) %>%
    arrange(Day) %>%
    mutate(
      cum_disease = cumsum(nr_disease),
      cum_dropout = cumsum(nr_dropout)
    ) %>%
    ungroup() %>%
    # compute nr_at_risk, nr_survived and the survival proportion within each interval
    mutate(
      nr_at_risk  = nr_enrolled - cum_dropout - (cum_disease - nr_disease),
      nr_survived = nr_enrolled - cum_dropout - cum_disease,
      # used to compute the Kaplan-Meier estimator
      interval_survival = nr_survived / nr_at_risk,
      # used to compute the variance of the Kaplan-Meier estimator using Greenwood's formula
      interval_variance = nr_disease / (nr_at_risk * nr_survived)
    ) %>%
    # compute KM estimator
    group_by(Location, Drug_regimen) %>%
    arrange(Day) %>%
    mutate(
      # Kaplan-Meier estimator
      cum_survival = cumprod(interval_survival),
      # variance of Kaplan-Meier estimator using Greenwood's formula
      var_survival = cum_survival^2 * cumsum(interval_variance)
    ) %>%
    ungroup() %>%
    # cleanup
    mutate(
      cum_survival = if_else(nr_survived > 0, cum_survival, NA),
      var_survival = if_else(nr_survived > 0, var_survival, NA)
    ) %>%
    # compute cumulative risk
    mutate(
      cum_risk = 100 * (1 - cum_survival),
      sd_risk  = 100 * sqrt(var_survival)
    )
  
  
  dt_observed <-
    dt_observed %>%
    # add Day == 0 for Location & Drug_regimen combinations that are without
    bind_rows(
      
      setdiff(
        expand_grid(
          Location     = unique(dt_observed$Location),
          Drug_regimen = unique(dt_observed$Drug_regimen),
          Day          = 0L
        ),
        dt_observed %>%
          filter(Day == 0L) %>%
          select(Location, Drug_regimen, Day)
      )
      
    ) %>%
    # fill missing value for cum_risk
    mutate(cum_risk = if_else(Day == 0 & is.na(cum_risk), 0, cum_risk)) %>%
    arrange(Location, Drug_regimen, Day)
  
  
  return(dt_observed)
  
}


compute_density <- function (dt) {
  
  dt_temp <- dt %>% compute_summary_stats()
  
  dt %>%
    handle_burnin(action = "remove") %>%
    # remove really large outliers and compute densities
    group_by(name) %>%
    filter(value < quantile(value, probs = 0.99)) %>%
    reframe(
      x = density(value, n = 2^8)$x,
      y = density(value, n = 2^8)$y
    ) %>%
    # add summaries
    left_join(dt_temp, by = "name") %>%
    # determine position w.r.t. 95% CrI
    mutate(
      position =
        case_when(
          x < q1 ~ "Lower 2.5%",
          x > q2 ~ "Upper 2.5%",
          TRUE   ~ "Middle 95%"
        )
    ) %>%
    # for each median value m on the x-axis and each facet, find the corresponding y value
    mutate(abs_diff = abs(x - m)) %>%
    group_by(name) %>%
    mutate(y_median = if_else(abs_diff == min(abs_diff), y, NA_real_)) %>%
    ungroup()
  
}


add_known_efficacies <- function (dt, step, taskname) {
  
  dt$eff__No_8AQs <- 0

  if (step != "A") {
    
    dt_temp <-
      read_csv(file.path(dirX_inp(step), paste0("input-", taskname, ".csv"))) %>%
      filter(Drug_regimen != "No_8AQs" & !is.na(eff)) %>%
      select(Drug_regimen, eff)
    
    for (i in 1:nrow(dt_temp)) {
      var <- paste0("eff__", dt_temp$Drug_regimen[i])
      dt[[var]] <- dt_temp$eff[i]
    }
    
    
  }
  
  return(dt)
  
}


transform_for_Kaplan_Meier <- function (dt) {
  
  # set seed for sampling
  set.seed(0)
  
  
  # randomly sampled iterations (keep chain column for inner_join() below)
  dt <- dt %>% slice_sample(n = n_sample)
  
  
  # get parameters of the same type in a column
  dt_rbite <-
    dt %>%
    select(chain, i, starts_with("rbite")) %>% 
    { 
      if (any(grepl("rbite__", colnames(.)))) {
        pivot_longer(
          .,
          cols = starts_with("rbite"),
          names_prefix = "rbite__",
          names_to = "Location",
          values_to = "rbite"
        )
      } else {
        . # do not delete this line!
      }
    }
  
  dt_rrela <-
    dt %>%
    select(chain, i, starts_with("rrela")) %>%
    { 
      if (any(grepl("rrela__", colnames(.)))) {
        pivot_longer(
          .,
          cols = starts_with("rrela"),
          names_prefix = "rrela__",
          names_to = "Location",
          values_to = "rrela"
        )
      } else {
        . # do not delete this line!
      }
    }
  
  dt_eff <-
    dt %>%
    select(chain, i, starts_with("eff")) %>%
    pivot_longer(
      cols = starts_with("eff"),
      names_prefix = "eff__",
      names_to = "Drug_regimen",
      values_to = "eff"
    )
  
  
  # collect parameters together
  inner_join(
    dt_rbite, dt_rrela,
    by = intersect(names(dt_rbite), names(dt_rrela))
  ) %>%
    inner_join(dt_eff, by = c("chain", "i"), relationship = "many-to-many") %>%
    select(-chain, -i)
  
}


add_input_variables <- function (dt, step, taskname) {
  
  dt_temp <-
    read_csv(file.path(dirX_inp(step), paste0("input-", taskname, ".csv"))) %>%
    select(any_of(c(
      "Drug_regimen", "Location", "Day",
      "proph_BS", "gamma_L", "h_tP", "AON", "p_meta"
    ))) %>%
    mutate(Day = max(Day)) %>%
    rename(maxday = Day) %>%
    unique()

  dt %>% right_join(dt_temp)
  
}

