
# Code adapted from: Griffin (2015) Plos Computational Biology
g <- function (kappa, t) {
  
  gamma <- beta(0.5, kappa + 0.5) / pi
  temp <- ( 1 + cos(2 * pi * t) ) / 2
  res <- (temp^kappa) / gamma
  
  return(res)
  
}


# implement equation 6 from Griffin2015a
## compute g at t = dur_high_seas / 2
## g reaches its maximum at t = 0, and, since g is symmetric, it cuts the line
## y = 1 at the end of the high season, i.e. after dur_high_seas / 2
g_minus_1 <- function (kappa, dur_high_seas) {
  
  return(g(kappa, dur_high_seas/2) - 1)
  
}


# find the power parameter, kappa, that yields the desired duration of the high
# transmission season, dur_high_seas
find_kappa <- function (dur_high_seas) {
  
  # first guesses for the root, kappa
  kappa0 <- 1
  kappa1 <- 1
  
  # we want kappa0 to yield g_minus_1 < 0
  while (g_minus_1(kappa1, dur_high_seas) > 0) {
    kappa0 <- kappa1
    kappa1 <- 2 * kappa1
  }
  
  # we want kappa1 to yield g_minus_1 > 0
  while (g_minus_1(kappa0, dur_high_seas) < 0) {
    kappa1 <- kappa0
    kappa0 <- kappa0 / 2
  }
  
  # the real root will be between kappa0 and kappa1
  kappa <-
    uniroot(
      f             = g_minus_1,
      interval      = c(kappa0, kappa1),
      dur_high_seas = dur_high_seas
    )$root
  
  return(kappa)
  
}


# empty input file table
get_empty_table <- function (file_type) {
  
  if (file_type == "par") {
    tribble(
      # ---------------------------------------------------------------------- #
      ~name,                ~description,
      # ---------------------------------------------------------------------- #
      "N_part",             "Number_of_participants",
      "EIR_equil",          "equilibrium_EIR",
      "P_occup",            "Prop_occupational_men",
      "risk_occup",         "Occupational_risk_factor",
      "start_time",         "Start_time_in_years",
      "end_time",           "End_time_in_years",
      "burnin_time",        "Burn_in_time_in_years",
      "age_mean",           "mean_population_age",
      "age_max",            "maximum_population_age",
      "P_preg",             "Proportion_women_pregnant_or_breastfeeding",
      "dur_preg",           "Duration_of_pregnancy_and_breastfeeding",
      "age_0",              "age-dependent_biting_parameter",
      "rho_age",            "age-dependent_biting_parameter",
      "sig_het",            "heterogeneity_exposure",
      "bb",                 "mosquito_to_human_transmission",
      "c_PCR",              "human_to_mosquito_transmission_PCR_detectable",
      "c_LM",               "human_to_mosquito_transmission_LM_detectable",
      "c_D",                "human_to_mosquito_transmission_clinical_disease",
      "c_T",                "human_to_mosquito_transmission_treatment",
      "d_latent",           "latent_duration_in_liver",
      "r_LM",               "LM_recovery_rate",
      "r_D",                "clinical_disease_recovery",
      "r_T",                "treatment_clearance_rate",
      "d_PCR_min",          "duration_of_PCR_infection_full_immunity",
      "d_PCR_max",          "duration_of_PCR_infection_no_immunity",
      "A_PCR_50pc",         "anti_parasite_immunity_scale_parameter",
      "K_PCR",              "anti_parasite_immunity_shape_parameter",
      "u_par",              "anti_parasite_immunity_scale_parameter",
      "r_par",              "anti_parasite_immunity_decay",
      "a_par_rel",          "increment_of_anti_parasite_immunity_following_a_relapse",
      "phi_LM_max",         "prob_BS_infection_no_immunity",
      "phi_LM_min",         "prob_BS_infection_full_immunity",
      "phi_LM_rel",         "following_relapse__factor_of_phi_LM",
      "A_LM_50pc",          "BS_immunity_scale_parameter",
      "K_LM",               "BS_immunity_shape_parameter",
      "u_clin",             "clinical_immunity_scale_parameter",
      "r_clin",             "clinical_immunity_decay",
      "a_clin_rel",         "increment_of_clinical_immunity_following_a_relapse",
      "phi_D_max",          "prob_clinical_episode_no_immunity",
      "phi_D_min",          "prob_clinical_episode_max_immunity",
      "phi_D_rel",          "following_relapse__factor_of_phi_D",
      "A_D_50pc",           "clinical_immunity_scale_parameter",
      "K_D",                "clinical_immunity_shape_parameter",
      "P_mat",              "proportion_maternal_immunity",
      "d_mat",              "maternal_immunity_decay",
      "ff",                 "relapse_rate",
      "gamma_L",            "liver_clearance_rate",
      "CM_regimen",         "treatment_regimen",
      "CM_cover",           "proportion_of_symptomatic_cases_treated",
      "CM_CQ_eff",          "chloroquine_treatment_efficacy",
      "CM_CQ_eff_wPQ",      "chloroquine_treatment_efficacy_withPQ",
      "CM_CQ_proph",        "chloroquine_treatment_prophylaxis",
      "CM_PQ_eff",          "primaquine_efficacy",
      "CM_PQ_proph",        "primaquine_prophylaxis",
      "CM_PQ_adhere",       "primaquine_adherence",
      "CM_PQ_lowage",       "minimum_age_primaquine",
      "CM_PQ_G6PD_risk",    "G6PD_risk_for_primaquine",
      "CM_PQ_CYP2D6_risk",  "CYP2D6_risk_for_primaquine",
      "CM_PQ_preg_risk",    "G6PD_risk_for_pregnancy",
      "CM_G6PD_test",       "G6PD_testing",
      "G6PD_prev",          "G6PD_deficiency_prevalence",
      "mu_G6PD_nor",        "mean_G6PD_activity_in_normals",
      "sig_G6PD_nor",       "standard_deviation_in_G6PD_activity_in_normals",
      "mu_G6PD_het",        "mean_G6PD_activity_in_heterozygous_deficient",
      "sig_G6PD_het",       "standard_deviation_in_G6PD_activity_in_heterozygous_deficient",
      "mu_G6PD_def",        "mean_G6PD_activity_in_homozygous_deficient",
      "sig_G6PD_def",       "standard_deviation_in_G6PD_activity_in_homozygous_deficient",
      "CYP2D6_prev",        "CYP2D6_prevalence",
      "rho_round_LLIN",     "correlation_between_LLIN_rounds",
      "rho_round_IRS",      "correlation_between_IRS_rounds",
      "rho_round_MDA",      "correlation_between_MDA_rounds",
      "rho_LLIN_IRS",       "correlation_between_IRS_LLINS",
      "rho_MDA_VC",         "correlation_between_MDA_vector_control"
    )
  } else if (file_type == "mosq") {
    tribble(
      # ----------------------------------------------------------------- #
      ~name,            ~description,
      # ----------------------------------------------------------------- #
      "mu_M",           "mosquito_death_rate",
      "tau_M",          "duration_of_sporogony",
      "d_E_larvae",     "early_larval_development",
      "d_L_larvae",     "late_larval_development",
      "d_pupae",        "pupal_development",
      "mu_E0",          "early_larval_death_rate",
      "mu_L0",          "late_larval_death_rate",
      "mu_P",           "pupal_death_rate",
      "beta_larvae",    "mosquito_egg_rate",
      "gamma_larvae",   "relative_density_dependence",
      "dry_seas",       "dry_season_proportion",
      "kappa_seas",     "seasonality_shape_parameter",
      "t_peak_seas",    "seasonal_peak_time",
      "Q_0",            "human_blood_index",
      "CHI_endo",       "endophily",
      "PSI_indoors",    "proportion_indoor_bites",
      "PSI_bed",        "proportion_bites_in_bed",
      "delta_1",        "foraging_time",
      "delta",          "gonotrophic_cycle",
      "LLIN_half_life", "LLIN_half_life",
      "PYR_half_life",  "pyrethroid_half_life",
      "r_LLIN_0",       "LLIN_insecticide_repellency",
      "r_LLIN_net",     "net_repellency",
      "d_LLIN_0",       "LLIN_death",
      "IRS_half_life",  "IRS_insecticide_half_life",
      "r_IRS_0",        "IRS_repellency",
      "d_IRS_0",        "IRS_death"
    )
  } else {
    tibble(
      name = c(
        "years",
        "CM0_cover",
        "CM0_eff",
        "CM0_proph",
        "CM1_cover",
        "CM1_CQ_eff",
        "CM1_CQ_eff_wPQ",
        "CM1_CQ_proph",
        "CM1_PQ_eff",
        "CM1_PQ_proph",
        "CM1_PQ_adhere",
        "CM1_PQ_lowage",
        "CM1_PQ_G6PD_risk",
        "CM1_PQ_CYP2D6_risk",
        "CM1_PQ_preg_risk",
        "CM1_G6PD_test",
        "CM2_cover",
        "CM2_CQ_eff",
        "CM2_CQ_eff_wPQ",
        "CM2_CQ_proph",
        "CM2_PQ_eff",
        "CM2_PQ_proph",
        "CM2_PQ_adhere",
        "CM2_PQ_lowage",
        "CM2_PQ_G6PD_risk",
        "CM2_PQ_CYP2D6_risk",
        "CM2_PQ_preg_risk",
        "CM2_TQ_eff",
        "CM2_TQ_proph",
        "CM2_TQ_adhere",
        "CM2_TQ_lowage",
        "CM2_TQ_G6PD_risk",
        "CM2_TQ_CYP2D6_risk",
        "CM2_TQ_preg_risk",
        "CM2_G6PD_test",
        "LLIN_cover",
        "IRS_cover",
        "MDA0_cover",
        "MDA0_CQ_eff",
        "MDA0_CQ_proph",
        "MDA1_cover",
        "MDA1_CQ_eff",
        "MDA1_CQ_eff_wPQ",
        "MDA1_CQ_proph",
        "MDA1_PQ_eff",
        "MDA1_PQ_proph",
        "MDA1_PQ_adhere",
        "MDA1_PQ_lowage",
        "MDA1_PQ_G6PD_risk",
        "MDA1_PQ_CYP2D6_risk",
        "MDA1_PQ_preg_risk",
        "MDA1_G6PD_test",
        "MDA2_cover",
        "MDA2_CQ_eff",
        "MDA2_CQ_eff_wPQ",
        "MDA2_CQ_proph",
        "MDA2_PQ_eff",
        "MDA2_PQ_proph",
        "MDA2_PQ_adhere",
        "MDA2_PQ_lowage",
        "MDA2_PQ_G6PD_risk",
        "MDA2_PQ_CYP2D6_risk",
        "MDA2_PQ_preg_risk",
        "MDA2_TQ_eff",
        "MDA2_TQ_proph",
        "MDA2_TQ_adhere",
        "MDA2_TQ_lowage",
        "MDA2_TQ_G6PD_risk",
        "MDA2_TQ_CYP2D6_risk",
        "MDA2_TQ_preg_risk",
        "MDA2_G6PD_test",
        "MSAT0_cover",
        "MSAT0_RDT_PCR",
        "MSAT0_sens",
        "MSAT0_CQ_eff",
        "MSAT0_CQ_proph",
        "MSAT1_cover",
        "MSAT1_RDT_PCR",
        "MSAT1_sens",
        "MSAT1_CQ_eff",
        "MSAT1_CQ_eff_wPQ",
        "MSAT1_CQ_proph",
        "MSAT1_PQ_eff",
        "MSAT1_PQ_proph",
        "MSAT1_PQ_adhere",
        "MSAT1_PQ_lowage",
        "MSAT1_PQ_G6PD_risk",
        "MSAT1_PQ_CYP2D6_risk",
        "MSAT1_PQ_preg_risk",
        "MSAT1_G6PD_test",
        "MSAT2_cover",
        "MSAT2_RDT_PCR",
        "MSAT2_sens",
        "MSAT2_CQ_eff",
        "MSAT2_CQ_eff_wPQ",
        "MSAT2_CQ_proph",
        "MSAT2_PQ_eff",
        "MSAT2_PQ_proph",
        "MSAT2_PQ_adhere",
        "MSAT2_PQ_lowage",
        "MSAT2_PQ_G6PD_risk",
        "MSAT2_PQ_CYP2D6_risk",
        "MSAT2_PQ_preg_risk",
        "MSAT2_TQ_eff",
        "MSAT2_TQ_proph",
        "MSAT2_TQ_adhere",
        "MSAT2_TQ_lowage",
        "MSAT2_TQ_G6PD_risk",
        "MSAT2_TQ_CYP2D6_risk",
        "MSAT2_TQ_preg_risk",
        "MSAT2_G6PD_test",
        "STAT1_cover",
        "STAT1_sens",
        "STAT1_spec",
        "STAT1_RDT_PCR",
        "STAT1_CQ_eff",
        "STAT1_CQ_eff_wPQ",
        "STAT1_CQ_proph",
        "STAT1_PQ_eff",
        "STAT1_PQ_proph",
        "STAT1_PQ_adhere",
        "STAT1_PQ_lowage",
        "STAT1_PQ_G6PD_risk",
        "STAT1_PQ_CYP2D6_risk",
        "STAT1_PQ_preg_risk",
        "STAT1_G6PD_test",
        "STAT2_cover",
        "STAT2_RDT_PCR",
        "STAT2_sens",
        "STAT2_spec",
        "STAT2_CQ_eff",
        "STAT2_CQ_eff_wPQ",
        "STAT2_CQ_proph",
        "STAT2_PQ_eff",
        "STAT2_PQ_proph",
        "STAT2_PQ_adhere",
        "STAT2_PQ_lowage",
        "STAT2_PQ_G6PD_risk",
        "STAT2_PQ_CYP2D6_risk",
        "STAT2_PQ_preg_risk",
        "STAT2_TQ_eff",
        "STAT2_TQ_proph",
        "STAT2_TQ_adhere",
        "STAT2_TQ_lowage",
        "STAT2_TQ_G6PD_risk",
        "STAT2_TQ_CYP2D6_risk",
        "STAT2_TQ_preg_risk",
        "STAT2_G6PD_test",
        "IVM_coverage",
        "IVM_baseline_efficacy",
        "IVM_half_live"
      )
    ) 
  }
  
}
