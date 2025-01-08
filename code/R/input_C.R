create_input_files_C <- function (step, CYP2D6_prev, proph_CQ, proph_TQl, proph_TQh, d_gamma) {
  
  # Input data -----------------------------------------------------------------
  # Data from a meta-analysis of clinical trials by 
  # [Watson et al.](https://doi.org/10.7554/eLife.83433)
  
  dt_trial <- 
    tribble(
      # ----------------------------------------------- #
      ~Drug_regimen, ~nr_enrolled, ~Day, ~nr_recurrences,
      # ----------------------------------------------- #
      "No_8AQs",              182,  120,             101,
      "PQ_lowdose14",         257,  120,              57,
      "TQ_lowdose",           368,  120,              79,
      "TQ_highdose",           54,  120,               4
    )
  
  # write to file
  write_csv(dt_trial, file = file.path(dirX_inp(step), "data.csv"))
  
  
  # Vary parameters ------------------------------------------------------------
  
  v0 <- 0.05
  v1 <- 14L
  v2 <- 100L
  
  dt_vary <- 
    tribble(
      # ------------------------------------------------------------------------------------- #
      ~taskname,     ~AON, ~CYP2D6_prev,   ~d_gamma,   ~proph_CQ,   ~proph_TQl,   ~proph_TQh,   
      # ------------------------------------------------------------------------------------- #
      "AON",         1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl,    proph_TQh,
      "leaky",       0,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl,    proph_TQh,
      "CYP2D6minus", 1,    CYP2D6_prev-v0, d_gamma,    proph_CQ,    proph_TQl,    proph_TQh,
      "CYP2D6plus",  1,    CYP2D6_prev+v0, d_gamma,    proph_CQ,    proph_TQl,    proph_TQh,
      "dgminus",     1,    CYP2D6_prev,    d_gamma-v2, proph_CQ,    proph_TQl,    proph_TQh,
      "dgplus",      1,    CYP2D6_prev,    d_gamma+v2, proph_CQ,    proph_TQl,    proph_TQh,
      "CQminus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ-v1, proph_TQl,    proph_TQh,
      "CQplus",      1,    CYP2D6_prev,    d_gamma,    proph_CQ+v1, proph_TQl,    proph_TQh,
      "TQlminus",    1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl-v1, proph_TQh,
      "TQlplus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl+v1, proph_TQh,
      "TQhminus",    1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl,    proph_TQh-v1,
      "TQhplus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_TQl,    proph_TQh+v1
    )
  
  # write to file
  write_csv(dt_vary, file = file.path(dirX_inp(step), "varying.csv"))
  
  
  # Loop over varying parameters -----------------------------------------------
  
  for (i in 1:nrow(dt_vary)) {
    
    row <- dt_vary[i, ]
    
    
    # Read median efficacy of PQ lowdose
    mytaskname <- if (grepl("TQ", row$taskname)) "AON" else row$taskname
    eff__PQ_lowdose <-
      read_csv(file.path(
        dirX_res("B"),
        paste0("MCMC-", mytaskname, ".csv")
      )) %>%
      # remove burnin
      filter(i > burnin * max(i)) %>%
      # select efficacy columns
      select(eff__PQ_lowdose) %>%
      summarise(eff__PQ_lowdose = median(eff__PQ_lowdose)) %>%
      pull()
    
    
    # Prep MCMC input data
    dt_trial %>%
      bind_cols(row) %>%
      mutate(
        # add BS prophylaxis period
        proph_BS =
          case_match(
            Drug_regimen,
            "TQ_lowdose"  ~ proph_TQl,
            "TQ_highdose" ~ proph_TQh,
            .default      = proph_CQ
          ),
        # add efficacy estimates
        eff =
          case_match(
            Drug_regimen,
            "No_8AQs"      ~ 0,
            "PQ_lowdose14" ~ eff__PQ_lowdose,
            .default       = NA
          ),
        # define hypnozoite clearance rate
        gamma_L = 1 / d_gamma,
        # remove BS prophylaxis from time variable
        delta_t = Day - proph_BS,
        # precompute h_tP
        h_tP = exp( -gamma_L * proph_BS),
        # probability of metabolising the liver-stage drug
        p_meta = if_else(grepl("PQ", Drug_regimen), 1 - CYP2D6_prev, 1)
      ) %>%
      # write to file
      write_csv(file = file.path(dirX_inp(step), paste0("input-", row$taskname, ".csv")))
    
  }
  
}


# Define names of parameter to be fit
# - rbite = daily reinfection (biting) rate
# - rrela = daily relapse rate
# - eff   = efficacy
define_parnames_C <- function () {
  
  drug_regimen <- c("TQ_lowdose", "TQ_highdose")
  par <- c("rbite__C", "rrela__C", "eff__TQ_lowdose", "eff__TQ_highdose")
  
  # get parameter ids
  id_rbite <- grep("rbite__", par)
  id_rrela <- grep("rrela__", par)
  id_eff   <- grep("eff__",  par)
  
  
  # return
  return(list(
    drug_regimen = drug_regimen,
    par          = par,
    id_rbite     = id_rbite,
    id_rrela     = id_rrela,
    id_eff       = id_eff
  ))
  
}

