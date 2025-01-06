fig_LHS_and_MCMC <- function (step, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": LHS and MCMC summary\n"))
  
  # fig name
  fig_name <- paste(step, "LHS_and_MCMC", sep = "-")
  
  
  # prep for plotting
  DT <- list()
  for (filetype in c("LHS", "MCMC")) {
    
    dt <-
      read_csv(file.path(dirX_inp(step), "varying.csv")) %>%
      group_split(taskname) %>%
      lapply(function (row) {
        
        # read in
        dt_temp <-
          read_estimates(step, filetype, row$taskname) %>% 
          pivot_and_rescale()
        
        if (filetype == "LHS") {
          dt_temp <- 
            dt_temp %>%
            filter(i %in% 1:n_chains) %>% 
            rename(m = value) %>%
            mutate(q1 = m, q2 = m)
        } else {
          dt_temp <- dt_temp %>% compute_summary_stats()
        }
        
        bind_cols(row, dt_temp, filetype = filetype)
        
      }) %>%
      bind_rows() %>%
      mutate(taskname = factor(taskname))
    
    DT <- c(DT, list(dt))
    
  }
  
  dt <- bind_rows(DT)
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = m, xmin = q1, xmax = q2, y = taskname)
    ) +
    geom_errorbar(colour = "grey60", width = ifelse(filetype == "LHS", 0, 1)) +
    geom_point() +
    facet_grid(rows = vars(filetype), cols = vars(name), scales = "free", switch = "y") +
    labs(x = "Task name", y = NULL) +
    theme_global(base_size = 8) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(size = 6.5),
      axis.text.y = element_text(hjust = 0),
      axis.ticks.y = element_blank()
    )
  
  
  # print to file
  my.ggsave(plot = p, path = dir_diag, filename = fig_name,
            device = device, width = length(unique(dt$name)) + 1, height = 3)
  
}


fig_MCMC_stats <- function (step, taskname, device = device_SI) {
  
  
  cat(paste0("# Fig ", step, ": MCMC stats, ", taskname, "\n"))
  
  # fig name
  fig_name <- paste(step, "MCMC_stats", taskname, sep = "-")
  
  p <- list()
  # acceptance rates --------------------------------------------------------- #
  
  # read in
  dt <-
    read_csv(file.path(dirX_res(step), paste0("MCMC_ar-", taskname, ".csv"))) %>%
    pivot_longer(!chain) %>%
    mutate(chain = as.factor(chain))
  
  
  # plot
  p[[1]] <-
    ggplot(
      data = dt,
      mapping = aes(x = chain, y = value, fill = chain)
    ) +
    geom_col() +
    facet_grid(cols = vars(name)) +
    ylim(0, NA) +
    labs(x = "Chain", y = "Acceptance rate (%)")
  
  
  # step size ---------------------------------------------------------------- #
  
  # read in
  dt <-
    read_csv(file.path(dirX_res(step), paste0("MCMC_steps-", taskname, ".csv"))) %>%
    thin_chains() %>%
    pivot_longer(!c(chain, i)) %>%
    mutate(chain = as.factor(chain))
  
  
  # plot
  p[[2]] <-
    ggplot(
      data = dt,
      mapping = aes(x = i, y = value, colour = chain)
    ) +
    geom_line() +
    facet_grid(cols = vars(name)) +
    scale_x_continuous(labels = number) +
    labs(x = "MCMC iteration", y = "Step size") +
    guides(colour = guide_legend(nrow = 1))
  
  
  # Autocorrelation ---------------------------------------------------------- #
  
  # read in
  dt <-
    read_estimates(step, "MCMC", taskname) %>%
    pivot_and_rescale() %>%
    # ignore log-likelihood
    filter(name != "ll") %>%
    # compute autocorrelation
    arrange(i) %>%
    group_by(name, chain) %>%
    reframe(lag  = acf (value, plot = FALSE)$lag[ , , 1],
            acf  = acf (value, plot = FALSE)$acf[ , , 1])
  
  
  # plot
  p[[3]] <-
    ggplot(
      data = dt,
      mapping = aes(x = lag, y = acf, colour = chain, fill = chain)
    ) +
    geom_col(width = 0.1, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed") +
    facet_grid(cols = vars(name)) +
    scale_y_continuous(name = "Autocorrelation", breaks = seq(-1, 1, by = 0.2))
  
  
  # combine plots --------------------------------------------------------------
  
  # add common theme
  p <- lapply(p, function (q) {q + theme_global(base_size = 8)})
  
  g <- my.arrange_panels(
    plotlist = p,
    add_tag = FALSE, 
    standardise_width = TRUE,
    common_legend = "bottom",
    common_legend_panel = 2
  )
  g <- arrangeGrob(grobs = g, heights = c(1,1,1, 0.2))
  
  
  # print to file
  my.ggsave(plot = g, path = dir_diag, filename = fig_name,
            device = device, width = length(unique(dt$name)) + 1, height = 5)
  
}


fig_MCMC_traces <- function (step, taskname, device = device_SI) {
  
  cat(paste0("# Fig ", step, ": MCMC traces, ", taskname, "\n"))
  
  # fig name
  fig_name <- paste(step, "traces", taskname, sep = "-")
  
  
  # read in
  dt <-
    read_estimates(step, "MCMC", taskname) %>%
    pivot_and_rescale()
  
  
  # transform and plot
  p <-
    list(
      dt %>%
        handle_burnin(action = "extract") %>%
        thin_chains() %>%
        mutate(chainpart = "Just burnin"),
      dt %>%
        handle_burnin(action = "remove") %>%
        thin_chains() %>%
        mutate(chainpart = "No burnin")
    ) %>%
    lapply(function (ddt) {
      
      ggplot(
        data = ddt,
        mapping = aes(x = i, y = value, colour = chain)
      ) +
        geom_line() +
        facet_wrap(vars(name), nrow = 1, scales = "free") +
        scale_x_continuous(name = "MCMC iteration", labels = number) +
        scale_y_continuous(name = unique(ddt$chainpart), labels = number) +
        guides(color = guide_legend(nrow = 1)) +
        theme_global(base_size = 8) +
        theme(
          strip.text = element_text(size = 5.5),
          axis.text.x = element_text(hjust = 1, angle = 45)
        )
      
    })
  
  
  # combine plots
  g <- my.arrange_panels(
    plotlist = p,
    add_tag = FALSE, 
    standardise_width = TRUE,
    common_legend = "bottom"
  )
  g <- arrangeGrob(grobs = g, heights = c(1,1, 0.2))
  
  
  # print to file
  my.ggsave(plot = g, path = dir_diag, filename = fig_name,
            device = device, width = length(unique(dt$name)) + 1, height = 4)
  
}


fig_data_A <- function (device = device_SI) {
  
  step <- "A"  
  cat(paste0("# Fig ", step, ": input data\n"))
  
  
  # fig name
  fig_name <- paste(step, "data", sep = "-")
  
  
  # read in
  dt <- read_data_A(step)
  
  
  # transform
  dt <-
    dt %>%
    # get missing days
    right_join(
      expand_grid(
        Location     = unique(dt$Location),
        Drug_regimen = unique(dt$Drug_regimen),
        Day          = min(dt$Day):max(dt$Day)
      )
    ) %>%
    # manually fill first time point
    group_by(Location, Drug_regimen) %>%
    mutate(
      is_min = if_else(Day == min(Day), TRUE, FALSE),
      nr_enrolled = max(nr_enrolled, na.rm = TRUE),
      cum_dropout = if_else(is_min & is.na(cum_dropout), 0, cum_dropout),
      cum_disease = if_else(is_min & is.na(cum_disease), 0, cum_disease),
      nr_survived = if_else(is_min, nr_enrolled - cum_dropout - cum_disease, nr_survived)
    ) %>%
    ungroup() %>%
    # pivot
    select(Drug_regimen, Location, Day, nr_survived, cum_disease, cum_dropout) %>%
    pivot_longer(cols = c(nr_survived, cum_disease, cum_dropout)) %>%
    # fill missing values
    group_by(Location, Drug_regimen, name) %>%
    arrange(Day) %>%
    fill(value) %>%
    ungroup()
  
  
  # plot
  p <-
    ggplot(
      data = dt,
      mapping = aes(x = Day, y = value, fill = name)
    ) +
    geom_col(position = "fill") +
    facet_grid(rows = vars(Drug_regimen), cols = vars(Location), labeller = global_labeller) +
    labs(y = "Proportion of enrolled patients") +
    scale_x_continuous(
      breaks = c(0, 90, 180, 270, 365),
      minor_breaks = NULL
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("green4", "orange2", "purple4"),
      breaks = c("cum_disease", "cum_dropout", "nr_survived"),
      labels = c("Cumulative proportion with Pv recurrence",
                 "Cumulative proportion who dropped out",
                 "Patients still in the trial")
    ) +
    theme_global(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(10, "points")
    )
  
  
  # save to file
  my.ggsave(plot = p, path = dir_diag, filename = fig_name,
            device = device, width = 7.5, height = 4)
  
}
