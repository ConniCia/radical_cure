create_input_files_B <- function (step, CYP2D6_prev, proph_CQ, d_gamma) {
  
  # Input data -----------------------------------------------------------------
  # Data from a meta-analysis of clinical trials by 
  # [Commons et al.](https://doi.org/10.1016/S1473-3099(23)00430-9)
  
  dt_trial <- 
    tribble(
      # --------------------------------------------------- #
      ~Drug_regimen, ~nr_enrolled, ~Day,    ~m,   ~q1,   ~q2,# ~nr_observed_first_recurrences,
      # --------------------------------------------------- #
      "No_8AQs",             1470,  180, 0.510, 0.482, 0.539,# 628,
      "PQ_lowdose",          2569,  180, 0.193, 0.169, 0.219,# 183,
      "PQ_highdose",         2811,  180, 0.081, 0.070, 0.094#, 234
    )
  
  # write to file
  write_csv(dt_trial, file = file.path(dirX_inp(step), "data.csv"))
  
  
  # Vary parameters ------------------------------------------------------------
  
  v0 <- 0.05
  v1 <- 14L
  v2 <- 100L
  
  dt_vary <- 
    tribble(
      # --------------------------------------------------------- #
      ~taskname,     ~AON, ~CYP2D6_prev,   ~d_gamma,   ~proph_CQ,   
      # --------------------------------------------------------- #
      "AON",         1,    CYP2D6_prev,    d_gamma,    proph_CQ,
      "leaky",       0,    CYP2D6_prev,    d_gamma,    proph_CQ,
      "CYP2D6minus", 1,    CYP2D6_prev-v0, d_gamma,    proph_CQ,
      "CYP2D6plus",  1,    CYP2D6_prev+v0, d_gamma,    proph_CQ,
      "dgminus",     1,    CYP2D6_prev,    d_gamma-v2, proph_CQ,
      "dgplus",      1,    CYP2D6_prev,    d_gamma+v2, proph_CQ,
      "CQminus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ-v1,
      "CQplus",      1,    CYP2D6_prev,    d_gamma,    proph_CQ+v1
    )
  
  # write to file
  write_csv(dt_vary, file = file.path(dirX_inp(step), "varying.csv"))
  
  
  # Loop over varying parameters -----------------------------------------------
  
  for (i in 1:nrow(dt_vary)) {
    
    row <- dt_vary[i, ]
    
    
    # Read median efficacy of PQ highdose
    eff__PQ_highdose <-
      read_csv(file.path(
        dirX_res("A"),
        paste0("MCMC-", row$taskname, ".csv")
      )) %>%
      # remove burnin
      filter(i > burnin * max(i)) %>%
      # select efficacy columns
      select(starts_with("eff")) %>%
      summarise(across(everything(), median)) %>%
      # 
      pivot_longer(cols = everything()) %>%
      summarise(eff__PQ_highdose = mean(value)) %>%
      pull()
    
    
    # Prep MCMC input data
    dt_trial %>%
      bind_cols(row) %>%
      mutate(
        proph_BS = proph_CQ,
        nr_recurrences = round(m * nr_enrolled),
        # define hypnozoite clearance rate
        gamma_L = 1 / d_gamma,
        # remove BS prophylaxis from time variable
        delta_t = Day - proph_BS,
        # precompute h_tP
        h_tP = exp( -gamma_L * proph_BS),
        # probability of metabolising the liver-stage drug
        p_meta = if_else(grepl("PQ", Drug_regimen), 1 - CYP2D6_prev, 1)
      ) %>%
      # add efficacy estimates
      mutate(
        eff =
          case_match(
            Drug_regimen,
            "No_8AQs"     ~ 0,
            "PQ_highdose" ~ eff__PQ_highdose,
            .default      = NA
          )
      ) %>%
      # write to file
      write_csv(file.path(dirX_inp(step), paste0("input-", row$taskname, ".csv")))
    
  }
  
}


# Define names of parameter to be fit
# - rbite = daily reinfection (biting) rate
# - rrela = daily relapse rate
# - eff   = efficacy
define_parnames_B <- function () {
  
  drug_regimen <- "PQ_lowdose"
  par <- c("rbite__B", "rrela__B", "eff__PQ_lowdose")
  
  # get parameter ids
  id_rbite <- grep("rbite__", par)
  id_rrela <- grep("rrela__", par)
  id_eff   <- grep("eff__",   par)
  
  
  # return
  return(list(
    drug_regimen = drug_regimen,
    par          = par,
    id_rbite     = id_rbite,
    id_rrela     = id_rrela,
    id_eff       = id_eff
  ))
  
}

