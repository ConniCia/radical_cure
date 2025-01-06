fit_BC <- function (step, input_file, N_LHS, N_MCMC) {
  
  # Prep -----------------------------------------------------------------------
  
  # read in
  dt <- read_csv(input_file)
  
  
  # define helper objects
  par_names <- if (step == "B") define_parnames_B() else define_parnames_C()
  taskname  <- unique(dt$taskname)
  ll_bigF   <- if (unique(dt$AON) > 0) AON_bigF else leaky_bigF
  
  
  # clean
  dt <- dt %>%
    # calculations are much faster using base::data.frame()
    as.data.frame() %>%
    # only keep columns used in LHS and MCMC
    select(Drug_regimen, nr_recurrences, nr_enrolled, delta_t, h_tP, gamma_L, eff, p_meta)
  
  
  # Run LHS BC ------------------------------------------------------------------
  
  cat("# ", taskname, "\n")
  cat("# LHS:", N_LHS, "iterations\n")
  
  
  if (N_LHS > 0) {
    
    # run LHS
    tm <-
      system.time(
        par_LHS <- LHS_BC(N_LHS, dt, par_names, ll_bigF)
      )
    cat("Elapsed time = ", round(tm[["elapsed"]] / 60, 1), " minutes\n\n")
    
    
    # save to file
    par_LHS %>%
      as_tibble() %>%
      slice_head(n = 100) %>%
      write_csv(file = file.path(dirX_res(step), paste0("LHS-", taskname, ".csv")))
    
  }
  
  
  # Run MCMC BC -----------------------------------------------------------------
  
  cat("# MCMC:", N_MCMC, "iterations\n")
  
  
  if (N_MCMC > 0) {
    
    # get results of LHS
    par_LHS <- read_csv(file.path(dirX_res(step), paste0("LHS-", taskname, ".csv")), n_max = n_chains)
    
    
    # mclapply() outputs a list whose elements are named like the elements of the
    # first mclapply() input argument. The mclapply() output will be fed to map():
    # in order for this function not to fail, its arguments need to be named lists
    id_chains <- setNames(1:n_chains, 1:n_chains)
    
    
    # run MCMC in parallel
    tm <-
      system.time(
        res <-
          mclapply(
            id_chains,
            MCMC_BC,
            N_MCMC, par_LHS, dt, par_names, ll_bigF,
            mc.cores = n_chains
          )
      )
    cat("Elapsed time = ", round(tm[["elapsed"]] / 60, 1), " minutes\n\n")
    
    
    # combine results from different chains
    lst <- transpose(res) %>% map(bind_rows)
    
    
    # print chains to file
    lst$m_MCMC %>%
      # thinning to save computer space
      filter(i %in% seq.int(0, N_MCMC, by = 10)) %>%
      write_csv(file = file.path(dirX_res(step), paste0("MCMC-", taskname, ".csv")))
    
    
    # print step sizes to file
    lst$m_step %>%
      write_csv(file = file.path(dirX_res(step), paste0("MCMC_steps-", taskname, ".csv")))
    
    
    # print acceptance rates to file
    lst$m_accepted %>%
      mutate(across(.col = !chain, .fns = ~ round(100 * .x / N_MCMC, 2))) %>%
      write_csv(file = file.path(dirX_res(step), paste0("MCMC_ar-", taskname, ".csv")))
    
  }
  
}


LHS_BC <- function (N_LHS, dt, par_names, ll_bigF) {
  
  # get nr of parameters to fit
  N_par <- length(par_names$par)
  
  
  # retrieve parameter ids
  id_rbite <- par_names$id_rbite
  id_rrela <- par_names$id_rrela
  id_eff   <- par_names$id_eff
  
  
  # retrieve parameter names
  names_par          <- par_names$par
  names_drug_regimen <- par_names$drug_regimen
  
  
  # define limits of parameter values
  par_min_rates <- 0
  par_max_rates <- 5 / 365
  par_min_eff   <- 0.5
  par_max_eff   <- 1
  
  
  # set seed for sampling
  set.seed(0)
  
  
  # generate values from unit hypercube
  par_LHS <- lhs::randomLHS(n = N_LHS, k = N_par, preserveDraw = TRUE)
  
  
  # transform unit hypercube to given limits
  par_LHS <-
    sapply(1:N_par, function(j) {
      if (j %in% c(id_rbite, id_rrela)) {
        qunif(p = par_LHS[, j], min = par_min_rates, max = par_max_rates)
      } else {
        qunif(p = par_LHS[, j], min = par_min_eff,   max = par_max_eff)
      }
    })
  
  
  # compute log-likelihood for each hypercube value
  ll <-
    mclapply(1:N_LHS, function(i) {
      ll_BC(par_LHS[i, ], dt, ll_bigF, id_rbite, id_rrela, id_eff, names_drug_regimen)
    }, mc.cores = n_chains) %>%
    simplify2array()
  
  
  # # inspect results
  # if (any(! is.finite(ll)))
  #   stop("LHS does not yield any numeric log-likelihood value!\n")
  
  
  # add ll
  par_LHS <- cbind(par_LHS, unname(ll))
  
  
  # add column names
  colnames(par_LHS) <- c(names_par, "ll")
  
  
  # sort by ll
  par_LHS <- as_tibble(par_LHS) %>% arrange(desc(ll))
  
  
  return(par_LHS)
  
}


MCMC_BC <- function (chain, N_MCMC, par_LHS, dt, par_names, ll_bigF) {
  
  # prepare --------------------------------------------------------------------
  
  N_step <- as.integer(burnin * N_MCMC / 2) # end iteration of step size adapter
  N_cov  <- as.integer(N_step / 2)          # end iteration of covariance matrix tuning
  n_cov  <- max(3, as.integer(N_cov / 10))  # start iteration of covariance matrix tuning
  
  
  # get nr of parameters to fit
  N_par <- length(par_names$par)
  
  
  # retrieve parameter ids
  id_rbite <- par_names$id_rbite
  id_rrela <- par_names$id_rrela
  id_eff   <- par_names$id_eff
  
  
  # retrieve parameter names
  names_par          <- par_names$par
  names_drug_regimen <- par_names$drug_regimen
  
  
  # define limits of parameter values
  par_min       <- 0
  par_max_rates <- 15 / 365
  par_max_eff   <- 1
  
  
  # initialise parameter values
  par <- par_LHS[chain, 1:N_par] %>% unlist()
  names(par) <- names_par # superfluous, but doesn't hurt
  
  
  # initialise covariance matrix
  par_sd <- par/2
  
  
  # define step sizes and limits
  par_step <- rep(1, N_par)
  par_step_min <- 0.01
  par_step_max <- 15
  
  
  # matrix to hold changing parameter values
  m_MCMC <- matrix(NA, nrow = N_MCMC + 1, ncol = N_par + 2)
  colnames(m_MCMC) <- c(names_par, "ll", "i")
  
  
  # matrix to hold step sizes
  m_step <- matrix(NA, nrow = N_step + 1, ncol = N_par + 1)
  colnames(m_step) <- c(names_par, "i")
  
  
  # matrix to hold the nr of accepted proposals
  m_accepted <- matrix(0, nrow = 1, ncol = N_par)
  colnames(m_accepted) <- names_par
  
  
  # first MCMC step
  ll <- ll_BC(par, dt, ll_bigF, id_rbite, id_rrela, id_eff, names_drug_regimen)
  
  
  # stop if error
  if (! is.finite(ll))
    stop("MCMC starting point does not yield a numeric log-likelihood value!\n")
  
  
  # save first step in matrix
  i <- 1
  m_MCMC[i, ] <- c(par, ll, i)
  m_step[i, ] <- c(par_step, i)
  
  
  # run MCMC -------------------------------------------------------------------
  
  # set seed for sampling
  set.seed(0)
  
  
  # loop
  for (i in 2:(N_MCMC + 1)) {
    
    for (j in 1:N_par) {
      
      # propose new value
      new_value <- rnorm(n = 1, mean = par[j], sd = par_step[j] * par_sd[j])
      
      
      # only continue if proposed value is admissible
      ## (we assume a uniform prior distributions for all parameters: this
      ## yields a constant value for the prior probability, thus computing this
      ## would not add any value to the inference)
      if ( is_admissible(j, new_value, id_eff,
                         par_min, par_max_eff, par_max_rates) ) {
        
        # compute new ll
        new_par    <- par
        new_par[j] <- new_value
        new_ll <- ll_BC(new_par, dt, ll_bigF, id_rbite, id_rrela, id_eff, names_drug_regimen)
        delta <- new_ll - ll
        
        # decide whether to accept new values
        if (delta > 0 || runif(n = 1) < exp(delta)) {
          
          # accept new values
          par <- new_par
          ll  <- new_ll
          
          # update acceptance rate
          m_accepted[j] <- m_accepted[j] + 1
          
        }
        
        # adapt step size
        if (i <= N_step)
          par_step[j] <- adapt_par_step(i, delta, N_step,
                                        par_step[j], par_step_min, par_step_max)
        
        # update covariance matrix
        if (between(i, n_cov, N_cov))
          par_sd[j] <- sd(m_MCMC[1:(i-1), j])
        
      }
    } # end loop over j
    
    
    # # periodically print loop index
    # if (i %% 10 == 0 & chain == 1)
    #   cat("i =", i,
    #       "\nar   = ", paste(round(m_accepted*100/i, 1), collapse = "\t"),
    #       "\nstep = ", paste(round(par_step, 3),         collapse = "\t"),
    #       "\nsd   = ", paste(round(par_sd, 5),           collapse = "\t"),
    #       "\n\n")
    
    
    # save output in matrix
    m_MCMC[i, ] <- c(par, ll, i)
    if (i <= N_step + 1)
      m_step[i, ] <- c(par_step, i)
    
  } # end loop over i
  
  
  # return
  return(list(
    m_MCMC     = as_tibble(m_MCMC)     %>% mutate(chain = chain, i = i - 1),
    m_step     = as_tibble(m_step)     %>% mutate(chain = chain, i = i - 1),
    m_accepted = as_tibble(m_accepted) %>% mutate(chain = chain)
  ))
  
}


# Compute log-likelihood
ll_BC <- function (par, dt, ll_bigF, id_rbite, id_rrela, id_eff, names_drug_regimen) {
  
  dt %>%
    mutate(eff = replace(eff, Drug_regimen %in% names_drug_regimen, par[id_eff])) %>%
    # compute log-likelihood
    mutate(
      ll = 
        ll_binom(
          ll_bigF(delta_t, h_tP, gamma_L, par[id_rbite], par[id_rrela], 1 - eff * p_meta),
          nr_recurrences, 
          nr_enrolled
        ),
      .keep = "none"
    ) %>%
    colSums()
  
}


# binomial log-likelihood
# equivalent to: dbinom(x = k, size = n, prob = p, log = TRUE)
ll_binom <- function (p, k, n) {
  
  # binomial log-likelihood formula
  ## p = prop of patients still without Pv recurrence
  ## k = nr patients without Pv recurrence
  ## n = nr enrolled
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(1 - p)
  
}