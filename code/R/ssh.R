monitor_jobs <- function (compact = TRUE, ask = FALSE) {
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  
  # Check state
  if (compact) {
    ssh::ssh_exec_wait(
      session = session_maestro,
      command = c("squeue -u cciavare -r -o'%8j %T' | sort | uniq -c")
    )
  } else {
    ssh::ssh_exec_wait(
      session = session_maestro,
      command = c('squeue -u cciavare -o"%15i %9P %8j %8u %2t %7M %7l %5D %4C" --sort="t,i"')
    )
  }
  
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
  
  # Ask whether to continue
  if (ask) {
    
    res <- readline("Continue with analysis? (y/n) ")
    
    if (res == "n") {
      stop("Analysis interrupted by user.\n\n")
    } else {
      cat("Proceeding with analysis.\n\n")
    }
    
  }
  
}


cancel_jobs <- function () {
  
  # Get vector of currently running cluster jobs
  current_jobs <- monitor_jobs() %>% 
    capture.output() %>%
    paste(collapse = " ") %>%
    stringr::str_extract_all("\\w*RP\\w*") %>%
    unlist()
  
  
  if (length(current_jobs) < 1) {
    cat("No jobs to cancel since none running.\n")
  } else {
    
    # Open connection
    session_maestro <- ssh::ssh_connect(address_cluster)
    
    # loop through jobs to cancel
    for (jobname in current_jobs) {
      
      res <- readline(paste0("Cancel ", jobname, "? (y/n) "))
      
      if (res == "y") {
        ssh::ssh_exec_wait(
          session = session_maestro,
          command = paste0("scancel -u cciavare -n ", jobname)
        )
      }
      
    }
    
    # Close connection
    ssh::ssh_disconnect(session_maestro)
    
  }
  
}


wipe_workspace <- function (step) {
  
  # clean up local files: ./<dirX>
  res <- readline(paste0("Really delete LOCAL files for step ", step, "? (y/n) "))
  
  if (res == "y") {
    system(paste0("rm -rf ", dirX(step)))
    create_dirs(steps)
  }
  
  
  # clean up remote files: <path_on_cluster>/<dirX>
  res <- readline(paste0("Really delete REMOTE files for step ", step, "? (y/n) "))
  
  if (res == "y") {
    
    # Open connection
    session_maestro <- ssh::ssh_connect(address_cluster)
    
    
    ssh::ssh_exec_wait(
      session = session_maestro,
      command = paste0("rm -rf ", path_on_cluster, dirX(step))
    )
    
    # Close connection
    ssh::ssh_disconnect(session_maestro)
  }
  
  
  # clean up slurm files?
  res <- readline(paste0("Clean up slurm files? (y/n) "))
  
  if (res == "y") {
    
    # Open connection
    session_maestro <- ssh::ssh_connect(address_cluster)
    
    # clean up remote slurm log files: <home_on_cluster>/slurm*
    ssh::ssh_exec_wait(
      session = session_maestro,
      command = "find ~ -name 'slurm*.out' -type f -delete"
    )
    
    # Close connection
    ssh::ssh_disconnect(session_maestro)
    
  }
  
}


copy_input_to_cluster <- function (step, replace_code = "ask", replace_input = "ask") {
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  # Create input folders on remote volume
  ssh::ssh_exec_wait(
    session = session_maestro,
    command = paste0("mkdir -p ", path_on_cluster, dirX_inp(step))
  )
  
  
  # Copy and replace code
  if (replace_code != "no") {
    
    res <- readline(paste0("Really replace REMOTE code files? (y/n) "))
    
    if (res == "y") {
      
      # Delete existing code files (if any)
      ssh::ssh_exec_wait(
        session = session_maestro,
        command = paste0("rm -rf ", path_on_cluster, "code")
      )
      
      # Copy code files to server
      ssh::scp_upload(
        session = session_maestro,
        files   = "code",
        to      = paste0(path_on_cluster, "code"),
        verbose = FALSE
      )
      
      # Make bash files executable
      ssh::ssh_exec_wait(
        session = session_maestro,
        command = paste0("chmod +x ", path_on_cluster, "code/bash/*.sh")
      )
      
    }
    
  }
  
  
  # Copy and replace input files
  if (replace_input != "no") {
    
    res <- readline(paste0("Really replace REMOTE input files? (y/n) "))
    
    if (res == "y") {
      
      # Delete existing input files (if any)
      ssh::ssh_exec_wait(
        session = session_maestro,
        command = paste0("rm -rf ", path_on_cluster, dirX_inp(step))
      )
      
      # Copy input files to server
      ssh::scp_upload(
        session = session_maestro,
        files   = dirX_inp(step),
        to      = paste0(path_on_cluster, dirX_inp(step)),
        verbose = FALSE
      )
      
    }
    
  }
  
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
}


run_fit <- function (step, N_LHS, N_MCMC) {
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  # Launch jobs
  ssh::ssh_exec_wait(
    session = session_maestro,
    command = c(
      paste0("chmod +x ", path_on_cluster, "code/bash/*.sh"),
      paste(paste0(path_on_cluster, "code/bash/run_fit.sh"),
            path_on_cluster, step, N_LHS, N_MCMC, n_chains, burnin,
            sep = " ")
    )
  )
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
}


run_IBM <- function (step) {
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  
  # Make bash files executable
  ssh::ssh_exec_wait(
    session = session_maestro,
    command = paste0("chmod +x ", path_on_cluster, "code/bash/*.sh")
  )
  
  
  # (Re)compile C++ executable?
  if (file.exists(file.path(path_on_cluster, dirX_inp(step), "RP13_executable"))) {
    res <- readline("Recompile C++ executable? (y/n) ")
  } else {
    res <- "y"
  }
  
  if (res == "y") {
    cat("Compiling C++ executable.\n")
    
    ssh::ssh_exec_wait(
      session = session_maestro,
      command = paste(paste0(path_on_cluster, "code/bash/compile_IBM.sh"),
                      path_on_cluster, step,
                      sep = " ")
    )
  }
  
  
  # Launch jobs
  ssh::ssh_exec_wait(
    session = session_maestro,
    command = paste(paste0(path_on_cluster, "code/bash/run_IBM.sh"),
                    path_on_cluster, step,
                    sep = " ")
  )
  
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
}


run_clean <- function (step) {
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  
  # Launch jobs
  ssh::ssh_exec_wait(
    session = session_maestro,
    command =  c(
      paste0("chmod +x ", path_on_cluster, "code/bash/*.sh"),
      paste(paste0(path_on_cluster, "code/bash/clean_IBM.sh"),
            path_on_cluster, step,
            sep = " ")
    )
  )
  
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
}


copy_results_to_local <- function (steps_to_copy) {
  
  
  # Open connection
  session_maestro <- ssh::ssh_connect(address_cluster)
  
  
  # Get vector of currently non-empty results directories
  steps_with_results <- 
    sapply(steps_to_copy, function (step_tc) {
      
      nr_files <-
        ssh::ssh_exec_wait(
          session = session_maestro,
          command = paste0("ls ", file.path(path_on_cluster, dirX_res(step_tc)), " | wc -l")
        ) %>%
        capture.output() %>%
        first() %>%
        as.numeric()
      
      if (nr_files < 1) NULL else step_tc
      
    }) %>%
    unlist() %>%
    unname()
  
  
  if (length(steps_with_results) < 1) {
    cat("No results to copy since none exist.\n")
  } else {
    
    # loop through jobs whose results to copy to local
    for (step_wr in steps_with_results) {
      
      res <- readline(paste0("Overwrite results from ", step_wr, "? (y/n) "))
      
      if (res == "y") {
        
        # Compress results on cluster
        ssh::ssh_exec_wait(
          session = session_maestro,
          command = c(
            paste0("cd ", path_on_cluster),
            paste(
              "tar -czvf",
              # path of archive
              file.path(dirX(step_wr), "/results.tar.gz"),
              # path of what to archive
              "-C", dirX_res(step_wr), "."
            )
          )
        )
        
        # Copy compressed files from cluster to local
        ssh::scp_download(
          session = session_maestro,
          files   = file.path(path_on_cluster, dirX(step_wr), "/results.tar.gz"),
          to      = dirX(step_wr),
          verbose = FALSE
        )
        
        # Extract cleaned results on local
        system(paste0("rm -rf ", dirX_res(step_wr)))
        system(paste0("mkdir -p ", dirX_res(step_wr)))
        system(paste0("tar -xzvf ", dirX(step_wr), "/results.tar.gz -C ", dirX_res(step_wr)))
        
        # Clean up on remote
        ssh::ssh_exec_wait(
          session = session_maestro,
          command = paste0("rm ", path_on_cluster, dirX(step_wr), "/results.tar.gz")
        )
        
      }
    }
  }
  
  
  # Close connection
  ssh::ssh_disconnect(session_maestro)
  
}

