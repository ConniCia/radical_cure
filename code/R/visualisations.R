# A, B, C: diagnostic plots for MCMC
create_diagnostic_plots <- function () {
  
  # Data from IMPROV trial
  fig_data_A()
  
  
  for (step in c("A", "B", "C")) {
    
    # collective LHS and MCMC plot
    fig_LHS_and_MCMC(step)
    
    
    for (taskname in c("AON", "leaky")) {

      # MCMC stats
      # - acceptance rates
      # - step sizes
      # - autocorr
      fig_MCMC_stats(step, taskname)


      # MCMC traces
      # - just burnin
      # - no burnin
      fig_MCMC_traces(step, taskname)

    }
    
  }
  
  
}


## A, B, C: paper figures
create_figures_ABC <- function () {
  
  
  ## main: Kaplan-Meier: A
  fig_Kaplan_Meier(step = "A", taskname = "AON", device = device_main)
  
  
  ## SI: Kaplan-Meier: B, C
  fig_Kaplan_Meier(step = "B", taskname = "AON")
  fig_Kaplan_Meier(step = "C", taskname = "AON")
  
  
  ## SI: posterior densities: A, B, C
  fig_posterior_density(step = "A", taskname = "AON")
  fig_posterior_density(step = "A", taskname = "leaky")
  fig_posterior_density(step = "B", taskname = "AON")
  fig_posterior_density(step = "C", taskname = "AON")
  
  
  ## SI: sensitivity: A + B + C
  fig_sensitivity_fits()
  
  
  ## main: trial and operational conditions: A + B + C
  fig_survival_curves(step = "ABC", taskname = "AON", device = device_main)
  
}


## D: paper figures
create_figures_D <- function () {
  
  
  dt_timeline <- read_in(step = "D", file_name = "timeline.csv")
  dt_stats    <- read_in(step = "D", file_name = "stats.csv")
  dt_effsize  <- read_in(step = "D", file_name = "effsize.csv")
  
  
  year_endpoint <- 5
  
  
  ## main: trends moderate seasonality
  fig_trends(dt_timeline, year_endpoint, filter_seas = "Moderate", device = device_main)
  
  
  ## main: effsize
  fig_effsize(step = "D", dt_effsize, year_endpoint, device = device_main)
  
  
  ## SI: seasonality
  fig_seasonality(dt_timeline)
  
  
  ## SI: trends no seasonality
  fig_trends(dt_timeline, year_endpoint, filter_seas = "None")
  
  
  ## SI: trends high seasonality
  fig_trends(dt_timeline, year_endpoint, filter_seas = "High")
  
  
  ## SI: averted cases
  fig_cases_averted(dt_effsize, year_endpoint)
  
  
  ## SI: treatment courses
  fig_treatment_courses(dt_stats)
  
  
  ## SI: relapses and reinfections
  fig_relapses_reinfections(dt_stats)
  
  
  ## SI: sensitivity
  fig_sensitivity_IBM(dt_effsize, dt_stats, year_endpoint)
  
}


## A, B, C, D: paper tables
create_tables <- function () {
  
  ## input data
  # tab_eligibility_criteria()
  tab_drug_regimen_overview()
  tab_trial_data_overview()
  tab_KM()
  tab_trial_data_A()
  tab_trial_data_B()
  tab_trial_data_C()
  tab_varying_recurrence_model()
  tab_varying_IBM()
  tab_posterior_draws_IBM()
  tab_fixed_IBM()
  
  ## results: proportion of IBM population eligible for 8-AQ regimens
  tab_proportion_eligible_for_8AQ()
  
  ## results: summary of posterior parameter fits: A + B + C
  tab_estimates_rrela_rbite(taskname = "AON")
  tab_estimates_efficacy_effectiveness(taskname = "AON")
  
  
}
