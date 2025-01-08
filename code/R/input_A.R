create_input_files_A <- function (step, CYP2D6_prev, proph_CQ, proph_DP, d_gamma) {
  
  # Input data -----------------------------------------------------------------
  # De-identified data on the time to first P vivax recurrence from the 
  # [IMPROV clinical trial](http://dx.doi.org/10.1016/S0140-6736(19)31285-1)
  
  # character columns: Location, Partner_drug, Drug_regimen
  # integer columns: fail, Day
  dt_trial <- read_csv(file.path(dir_raw, "placeholder_IMPROV_data.csv"))
  
  
  # write to file
  write_csv(dt_trial, file = file.path(dirX_inp(step), "data.csv"))

  
  # Vary parameters ------------------------------------------------------------
  
  v0 <- 0.05
  v1 <- 14L
  v2 <- 100L
  
  dt_vary <- tribble(
    # -------------------------------------------------------------------------------- #
    ~taskname,     ~AON, ~CYP2D6_prev,   ~d_gamma,   ~proph_CQ,   ~proph_DP,   ~country,
    # -------------------------------------------------------------------------------- #
    "AON",         1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "",
    "leaky",       0,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "",
    "CYP2D6minus", 1,    CYP2D6_prev-v0, d_gamma,    proph_CQ,    proph_DP,    "",
    "CYP2D6plus",  1,    CYP2D6_prev+v0, d_gamma,    proph_CQ,    proph_DP,    "",
    "dgminus",     1,    CYP2D6_prev,    d_gamma-v2, proph_CQ,    proph_DP,    "",
    "dgplus",      1,    CYP2D6_prev,    d_gamma+v2, proph_CQ,    proph_DP,    "",
    "CQminus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ-v1, proph_DP,    "",
    "CQplus",      1,    CYP2D6_prev,    d_gamma,    proph_CQ+v1, proph_DP,    "",
    "DPminus",     1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP-v1, "",
    "DPplus",      1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP+v1, "",
    "Afghanistan", 1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "AF",
    "Ethiopia",    1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "ET",
    "Indonesia",   1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "ID",
    "Vietnam",     1,    CYP2D6_prev,    d_gamma,    proph_CQ,    proph_DP,    "VN"
  )
  
  # write to file
  write_csv(dt_vary, file = file.path(dirX_inp(step), "varying.csv"))
  
  
  # Prep MCMC input data -------------------------------------------------------
  
  for (i in 1:nrow(dt_vary)) {
    
    row <- dt_vary[i, ]
    
    dt_trial %>%
      # make more compact
      group_by(Location, Drug_regimen, Day, fail, Partner_drug) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      # add varying parameters
      bind_cols(row) %>%
      filter({row$country == "" | grepl(row$country, Location)}) %>%
      mutate(
        # define hypnozoite clearance rate
        gamma_L = 1 / d_gamma,
        # add BS prophylaxis period
        proph_BS = if_else(Partner_drug == "CQ", proph_CQ, proph_DP),
        # remove BS prophylaxis from time variable
        delta_t = Day - proph_BS,
        # precompute h_tP
        h_tP = exp( -gamma_L * proph_BS),
        # probability of metabolising the liver-stage drug
        p_meta = if_else(grepl("PQ", Drug_regimen), 1 - CYP2D6_prev, 1)
      ) %>%
      # exclude negative times
      filter(delta_t > 0) %>%
      # write to file
      write_csv(file.path(dirX_inp(step), paste0("input-", row$taskname, ".csv")))
    
  }
  
}


# Define names of parameter to be fit
# - rbite = daily reinfection (biting) rate, varies by location
# - rrela = daily relapse rate, varies by location
# - eff   = efficacy, varies by drug regimen
define_parnames_A <- function (dt) {
  
  location     <- unique(dt$Location)
  drug_regimen <- c("PQ_highdose7", "PQ_highdose14")
  
  par <-
    c(
      paste0("rbite__", location),
      paste0("rrela__", location),
      paste0("eff__",   drug_regimen)
    )
  
  # get parameter ids
  id_rbite <- grep("rbite__", par)
  id_rrela <- grep("rrela__", par)
  id_eff   <- grep("eff__",   par)
  
  
  # return
  return(list(
    location     = location,
    drug_regimen = drug_regimen,
    par          = par,
    id_rbite     = id_rbite,
    id_rrela     = id_rrela,
    id_eff       = id_eff
  ))
  
}
