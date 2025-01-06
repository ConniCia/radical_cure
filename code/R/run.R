run <- function (step, N_LHS = 0, N_MCMC = 0) {
  
  
  # check job status: ask if continue
  monitor_jobs(ask = TRUE)
  
  
  # wipe from local AND cluster
  wipe_workspace(step)
  
  
  # create input files
  cat(paste0("Step ", step, ": create input files\n"))
  if (step == "A") {
    
    create_input_files_A(
      step        = "A",
      CYP2D6_prev = dt_default$CYP2D6_prev,
      proph_CQ    = dt_default$CQ__proph,
      proph_DP    = dt_default$DP__proph,
      d_gamma     = dt_default$d_gamma
    )
    
  } else if (step == "B") {
    
    copy_results_to_local(steps_to_copy = "A")
    create_input_files_B(
      step        = "B",
      CYP2D6_prev = dt_default$CYP2D6_prev,
      proph_CQ    = dt_default$CQ__proph,
      d_gamma     = dt_default$d_gamma
    )
    
  } else if (step == "C") {
    
    copy_results_to_local(steps_to_copy = "B")
    create_input_files_C(
      step        = "C",
      CYP2D6_prev = dt_default$CYP2D6_prev,
      proph_CQ    = dt_default$CQ__proph,
      proph_TQl   = dt_default$TQ_lowdose__proph,
      proph_TQh   = dt_default$TQ_highdose__proph,
      d_gamma     = dt_default$d_gamma
    )
    
  } else if (step == "D") {
    
    copy_results_to_local(steps_to_copy = c("A", "B", "C"))
    create_input_files_D(step = "D")
    
  }
  
  
  # copy to cluster and make bash files executable
  cat(paste0("Step ", step, ": copy input files to cluster\n"))
  copy_input_to_cluster(step) 
  
  
  # run analysis on cluster
  cat(paste0("Step ", step, ": submit cluster jobs\n"))
  if (step %in% c("A", "B", "C")) {
    
    run_fit(step = step, N_LHS = N_LHS, N_MCMC = N_MCMC)
    
  } else if (step == "D") {
    
    run_IBM(step = "D")
    
  }
  
  
  # Monitor jobs
  monitor_jobs()
  
}


run_cleaning <- function (step) {
  
  
  # check job status: ask if continue
  monitor_jobs(ask = TRUE)
  
  # copy code to cluster and make executable
  copy_input_to_cluster(step, replace_input = "no") 
  
  # clean output from cluster analysis
  run_clean(step = "D")
 
  # Monitor jobs
  monitor_jobs()
  
}


copy_cluster_results <- function () {
  
  
  # check job status: ask if continue
  monitor_jobs(ask = TRUE)
  
  
  # copy results from cluster to local
  copy_results_to_local(steps_to_copy = steps)
  
}
