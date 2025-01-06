default_parameters <- function (print_to_file = FALSE) {
  
  dt <- tribble(
    # ------------------------------------------------------------------------ #
    ~Parameter,      ~Value, ~References, ~Description,
    # ------------------------------------------------------------------------ #
    
    # Default values of varying parameters -------------------------------------
    ##~Parameter,    ~Value, ~References, ~Description,
    
    "yEIR",          1,      "",          "Yearly EIR at equilibrium",
    "d_ff",          41,     "",          "Relapse periodicity (in days)",
    "dry_seas",      0.1,    "",          "Transmission level during dry season as compared to yearly average",
    "dur_high_seas", 0.5,    "",          "Proportion of the year that is dry season",
    "CM_cover",      0.9,    "",          "Proportion of symptomatic cases treated",
    "G6PD_prev",     0.05,   "",          "Prevalence of hemizygous G6PD deficient males (genotype)",
    "CYP2D6_prev",   0.05,   "",          "Prevalence of CYP2D6 low metabolisers",
    "phi_LM_rel",    1,      "",          "In case of a relapse - Factor for probability of LM-detectable infection",
    "phi_D_rel",     1,      "",          "In case of a relapse - Factor for probability of symptomatic disease",
    "a_par_rel",     1,      "",          "In case of a relapse - Factor for increment of anti-parasite immunity",
    "a_clin_rel",    1,      "",          "In case of a relapse - Factor for increment of clinical immunity",
    "risk_occup",    0,      "",          "Relative risk of occupational w.r.t. domestic exposure to mosquito",
    "P_occup",       0.5,    "",          "Proportion of adult men with additional mosquito exposure from occupation",
    
    
    # Parameters sampled from posterior distribution ---------------------------
    ##~Parameter,    ~Value, ~References,  ~Description,
    
    "sig_het",       NA,     "@White2018", "Heterogeneity in exposure, standard deviation on a log scale",
    "d_LM",          NA,     "@White2018", "Duration of LM-detectable infection",
    "d_PCR_max",     NA,     "@White2018", "Duration of PCR-detectable infection (no immunity)",
    "A_PCR_50pc",    NA,     "@White2018", "Anti-parasite immunity - scale parameter",
    "K_PCR",         NA,     "@White2018", "Anti-parasite immunity - shape parameter",
    "u_par",         NA,     "@White2018", "Acquisition of anti-parasite immunity - scale parameter",
    "phi_LM_max",    NA,     "@White2018", "Probability of LM-detectable infection (no immunity)",
    "phi_LM_min",    NA,     "@White2018", "Probability of LM-detectable infection (full immunity)",
    "A_LM_50pc",     NA,     "@White2018", "Anti-parasite immunity - scale parameter",
    "K_LM",          NA,     "@White2018", "Anti-parasite immunity - shape parameter",
    "u_clin",        NA,     "@White2018", "Acquisition of clinical immunity - scale parameter",
    "phi_D_max",     NA,     "@White2018", "Probability of clinical disease (no immunity)",
    "phi_D_min",     NA,     "@White2018", "Probability of clinical disease (full immunity)",
    "A_D_50pc",      NA,     "@White2018", "Acquisition of clinical immunity - scale parameter",
    "K_D",           NA,     "@White2018", "Acquisition of clinical immunity - shape parameter",
    "P_mat",         NA,     "@White2018", "New-born immunity relative to mother's",
    "d_mat",         NA,     "@White2018", "Duration of maternal immunity (in days)",
    
    
    # General parameters -------------------------------------------------------
    ##~Parameter,        ~Value,         ~References,                 ~Description,
    
    # simulation parameters (in years)
    "start_time",        -5,             "",                          "Simulation start time (in years)",
    "intervention_time", 0,              "",                          "Intervention start time (in years)",
    "end_time",          10.01,          "",                          "Simulation end time (in years)",
    "burnin_time",       10,             "",                          "Burn in time (in years)",
    
    # population size and demography
    "N_part",            20000,          "",                            "Number of participants",
    "age_mean",          22.5 * 365,     "Assumption",                  "Mean population age (in days)",
    "age_max",           80 * 365,       "Assumption",                  "Maximum population age (in days)",
    "prop_females",      0.5,            "Assumption",                  "Proportion of females in the population",
    "P_preg",            0.125,          "Assumption",                  "Proportion of adult women that are pregnant or breastfeeding a baby of <6 months",
    "age_preg_min",      16 * 365,       "Assumption",                  "Minimum age for start of pregnancy (in days)",
    "age_preg_max",      45 * 365,       "Assumption",                  "Maximum age for start of pregnancy (in days)",
    "dur_preg",          (9 + 6) * 30,   "Assumption",                  "Duration of pregnancy and breastfeeding (in days)",
    "age_0",             8 * 365,        "[@Carnevale1978; @Port1980]", "Age-dependent biting parameter",
    "rho_age",           0.85,           "[@Carnevale1978; @Port1980]", "Age-dependent biting parameter",
    
    # solve equilibrium ODE for baseline case management intervention: blood stage drug only (neither PQ nor TQ)
    "CM_regimen",        0,              "",                          "Case management treatment regimen (0 = blood stage drug only)",
    "CM_PQ_eff",         0,              "",                          "",
    "CM_PQ_proph",       0,              "",                          "",
    "CM_PQ_adhere",      0,              "",                          "",
    "CM_PQ_lowage",      0,              "",                          "",
    "CM_PQ_G6PD_risk",   0,              "",                          "",
    "CM_PQ_CYP2D6_risk", 0,              "",                          "",
    "CM_PQ_preg_risk",   0,              "",                          "",
    "CM_G6PD_test",      0,              "",                          "Test for G6PD deficiency before giving baseline treatment",
    
    # transmission dynamics (human and mosquito part)
    "bb",                0.5,            "Assumption",              "Probability of mosquito-to-human transmission",
    "c_PCR",             0.035,          "@Kiattibutr2017",         "Probability of (PCR-detectable) human-to-mosquito transmission",
    "c_LM",              0.10,           "@Kiattibutr2017",         "Probability of (LM-detectable) human-to-mosquito transmission",
    "c_D",               0.80,           "@Kiattibutr2017",         "Probability of human-to-mosquito transmission",
    "c_T",               0.40,           "@Kiattibutr2017",         "Probability of human-to-mosquito transmission",
    "d_latent",          10,             "@Herrera2011",            "Duration of latency in the liver (in days)",
    "r_D",               0.20,           "Assumption [@White2018]", "Rate of recovery from symptomatic disease",
    "r_T",               1.0,            "@Pukrittayakamee2008",    "Rate of progression through treatment",
    "d_PCR_min",         10,             "Assumption [@White2018]", "Duration of PCR-detectable infection (full immunity)",
    "r_par",             1 / (10 * 365), "Assumption [@White2018]", "Rate of decay of anti-parasite immunity",
    "r_clin",            1 / (30 * 365), "Assumption [@White2018]", "Rate of decay of clinical immunity",
    "d_gamma",           383,            "@White2014a",             "Mean duration of hypnozoite carriage (in days)",
    
    # parameters describing distribution of G6PD activity (phenotype) levels in population
    "mu_G6PD_nor",       10.16,          "@Nekkab2021",    "Mean G6PD activity - homozygous normal females and hemizygous normal males",
    "sig_G6PD_nor",      2.39,           "@Nekkab2021",    "Standard deviation in G6PD activity - homozygous normal females and hemizygous normal males",
    "mu_G6PD_het",       4.6,            "@Nekkab2021",    "Mean G6PD activity - heterozygous deficient females",
    "sig_G6PD_het",      1.67,           "@Nekkab2021",    "Standard deviation in G6PD activity - heterozygous deficient females",
    "mu_G6PD_def",       0.48,           "@Nekkab2021",    "Mean G6PD activity - homozygous deficient females and hemizygous deficient males",
    "sig_G6PD_def",      0.31,           "@Nekkab2021",    "Standard deviation in G6PD activity - homozygous deficient females and hemizygous deficient males",
    
    # parameters describing the sensitivity and specificity of SD Biosensor STANDARD G6PD test
    "G6PDtest_sens30",   100,            "@Pal2019",       "SD Biosensor STANDARD G6PD test - sensitivity in those with < 30% G6PD activity",
    "G6PDtest_spec30",   97,             "@Pal2019",       "SD Biosensor STANDARD G6PD test - specificity in those with < 30% G6PD activity",
    "G6PDtest_sens70",   95.5,           "@Pal2019",       "SD Biosensor STANDARD G6PD test - sensitivity in those with < 70% G6PD activity",
    "G6PDtest_spec70",   97,             "@Pal2019",       "SD Biosensor STANDARD G6PD test - specificity in those with < 70% G6PD activity",
    
    # unused parameters describing vector control or pop-level interventions:
    # correlation of pop samples who get successive rounds of interventions
    "rho_round_LLIN",    0.50,           "",            "Correlation between successive rounds of LLIN",
    "rho_round_IRS",     0.80,           "",            "Correlation between successive rounds of IRS",
    "rho_round_MDA",     0.50,           "",            "Correlation between successive rounds of MDA",
    "rho_LLIN_IRS",      0.60,           "",            "Correlation between successive rounds of LLIN and IRS",
    "rho_MDA_VC",        0.20,           "",            "Correlation between successive rounds of MDA and vector control",
    
    
    # Drug regimen parameters --------------------------------------------------
    ##~Parameter,            ~Value, ~References,           ~Description,
    
    # Blood stage drugs
    "CQ__eff",               1,      "Assumption",          "CQ efficacy",
    "CQ__eff_wPQ",           1,      "Assumption",          "CQ efficacy when given with PQ",
    "CQ__proph",             28,     "@Karunajeewa2008",    "Duration of CQ prophylaxis (in days)",
    "DP__proph",             42,     "@Karunajeewa2008",    "Duration of DP prophylaxis (in days)",
    
    # AQ ideal
    "AQ_ideal__eff",         1,      "Assumption",          "Magic bullet - hypnozoiticidal efficacy",
    "AQ_ideal__proph",       28,     "Assumption",          "Magic bullet - prophylaxis (in days)",
    "AQ_ideal__adhere",      1,      "Assumption",          "Magic bullet - probability of full adherence",
    "AQ_ideal__lowage",      0,      "Assumption",          "Magic bullet - lower age boundary for treatment prescription",
    "AQ_ideal__G6PD_risk",   0,      "Assumption",          "Magic bullet - presence of risk for G6PD deficient patients",
    "AQ_ideal__CYP2D6_risk", 0,      "Assumption",          "Magic bullet - ineffective for CYP2D6 low-metaboliser patients",
    "AQ_ideal__preg_risk",   0,      "Assumption",          "Magic bullet - presence of risk for pregnant and breastfeeding women",
    "AQ_ideal__G6PD_test",   0,      "Assumption",          "Magic bullet - test for G6PD deficiency before administering treatment",
    
    # PQ
    "PQ_7__proph",           8,      "Assumption",                 "PQ over 7 days - duration of prophylaxis (in days)",
    "PQ_14__proph",          8+7,    "Assumption",                 "PQ over 14 days - duration of prophylaxis (in days)",
    "PQ_7__adhere",          0.67,   "@Almeida2014",               "PQ over 7 days - adherence to full treatment course",
    "PQ_14__adhere",         0.57,   "Assumption",                 "PQ over 14 days - adherence to full treatment course",
    "PQ__lowage",            6*30,   "@WHO2023malariaguidelines",  "PQ - lower age boundary for treatment prescription (in days)",
    "PQ__G6PD_risk",         1,      "@WHO2023malariaguidelines",  "PQ - presence of risk for G6PD deficient patients",
    "PQ__CYP2D6_risk",       1,      "[@Bennett2013; @StJean2016; @Silvino2016]", "PQ - ineffective for CYP2D6 low-metaboliser patients",
    "PQ__preg_risk",         1,      "@WHO2023malariaguidelines",  "PQ - presence of risk for pregnant and breastfeeding women (baby's G6PD status unknown)",
    "PQ__G6PD_test",         1,      "@WHO2018G6PD",               "PQ - test for G6PD deficiency before administering treatment",
    
    # TQ
    "TQ_lowdose__proph",     45,     "@Lacerda2019",        "Low-dose TQ - duration of prophylaxis (in days)",
    "TQ_highdose__proph",    60,     "Assumption",          "High-dose TQ - duration of prophylaxis (in days)",
    "TQ__adhere",            1,      "Assumption",          "TQ - adherence to full treatment course",
    "TQ__lowage",            2*365,  "Assumption",          "TQ - lower age boundary for treatment prescription (in days)",
    "TQ__G6PD_risk",         1,      "@Rueangweerayut2017", "TQ - presence of risk for G6PD deficient patients",
    "TQ__CYP2D6_risk",       0,      "@StJean2016",         "TQ - ineffective for CYP2D6 low-metaboliser patients",
    "TQ__preg_risk",         1,      "@Rueangweerayut2017", "TQ - presence of risk for pregnant and breastfeeding women (baby's G6PD status unknown)",
    "TQ__G6PD_test",         1,      "@WHO2018G6PD",        "TQ - test for G6PD deficiency before administering treatment",
    
    
    # Mosquito parameters ------------------------------------------------------
    ##~Parameter,     ~Value, ~References,    ~Description,
    
    # mosquito population : An. gambiae s.l., Nigeria (Garki project)
    "mu_M",           1 / 6,  "@Nekkab2021",  "Mosquito death rate (daily)",
    "tau_M",          8.0,    "@Nekkab2021",  "Duration of sporogony (in days)",
    "d_E_larvae",     6.64,   "@White2011a",  "Duration of early larval instar stage (in days)",
    "d_L_larvae",     3.72,   "@White2011a",  "Duration of late larval instar stage (in days)",
    "d_pupae",        0.64,   "@White2011a",  "Duration of pupae stage (in days)",
    "mu_E0",          0.034,  "@White2011a",  "Mortality rate of early larval instars (when density is low)",
    "mu_L0",          0.035,  "@White2011a",  "Mortality rate of late larval instars (when density is low)",
    "mu_P",           0.25,   "@White2011a",  "Mortality rate of pupae",
    "beta_larvae",    21.19,  "@White2011a",  "Number of eggs laid per day per mosquito",
    "gamma_larvae",   13.25,  "@White2011a",  "Effect of density dependence on late instars w.r.t. early instars",
    "t_peak_seas",    0.0,    "",             "Timing of seasonal peak transmission period w.r.t. start of numerical simulation",
    
    # mosquito biting behaviour : An. darlingi, most common vector in Brazil
    "Q_0",            0.5,    "@Nekkab2021",  "Human blood index (proportion of blood meals taken on humans)",
    "CHI_endo",       0.9,    "@Nekkab2021",  "Proportion of endophilic mosquitoes (resting indoors after feeding)",
    "PSI_indoors",    0.2,    "@Nekkab2021",  "Proportion of bites taken on humans indoors",
    "PSI_bed",        0.2,    "@Nekkab2021",  "Proportion of bites taken on humans in bed",
    "delta_1",        0.68,   "@Nekkab2021",  "Time spent foraging for a blood meal",
    "delta",          3.0,    "@Killeen2000", "Time spent digesting blood meal (duration of gonotrophic cycle)",
    
    # effect of vector control on mosquitoes
    "LLIN_half_life", 1095.0, "",             "Half-life of loss of LLINs", 
    "PYR_half_life",  912.5,  "",             "Half-life of pyrethroid decay on LLINs", 
    "r_LLIN_0",       0.6,    "",             "Probability mosquito repelled (full insecticide activity)", 
    "r_LLIN_net",     0.2,    "",             "Probability mosquito repelled due to barrier effect of net (independent of insecticide)", 
    "d_LLIN_0",       0.3,    "",             "Probability mosquito dies during feeding attempt", 
    "IRS_half_life",  182.5,  "",             "Half-life of IRS insecticide decay", 
    "r_IRS_0",        0.5,    "",             "Probability mosquito repelled by IRS (full insecticide activity)", 
    "d_IRS_0",        0.4,    "",             "Probability mosquito killed by IRS (full insecticide activity)"
    
  )
  
  
  # return
  if (print_to_file) {
    dt
  } else {
    dt %>% filter(!is.na(Value)) %>% select(Parameter, Value) %>% pivot_wider(names_from = Parameter, values_from = Value)
  }
  
}


varying_parameters <- function () {
  
  bind_rows(
    
    # vary key parameters ---------------------------------------------------- #
    
    expand_grid(
      # yearly EIR at equilibrium
      yEIR = c(0.1, 1, 10),
      # adherence to PQ treatment courses
      tibble(
        PQ_7__adhere  = c(0.2, 0.67, 0.9),
        PQ_14__adhere = c(0.1, 0.57, 0.8)
      )
    ),
    
    expand_grid(
      # yearly EIR at equilibrium
      yEIR = c(0.1, 1, 10),
      # seasonality
      tibble(
        # EIR in dry season as proportion of mean yearly EIR
        dry_seas      = c(1  , 0.1, 0.05),
        # proportion of year with high transmission
        dur_high_seas = c(0.6, 0.5, 0.25)
      )
    ),
    
    # vary additional parameters --------------------------------------------- #
    
    # average duration between successive relapses
    tibble(d_ff = c(41, 69, 120)),
    # coverage of case management
    tibble(CM_cover = c(0.5, 0.7, 0.9, 1)),
    # prevalence of hemizygous G6PD deficient males
    tibble(G6PD_prev = c(0, 0.05, 0.1)),
    # prevalence of CYP2D6
    tibble(CYP2D6_prev = c(0, 0.05, 0.1)),
    # relative probability that relapse becomes LM-detectable w.r.t. reinfection
    tibble(phi_LM_rel = c(0.25, 0.5, 0.75, 1)),
    # relative probability that relapse becomes symptomatic w.r.t. reinfection
    tibble(phi_D_rel  = c(0.25, 0.5, 0.75, 1)),
    # relative increment of anti-parasite immunity following relapse w.r.t. reinfection
    tibble(a_par_rel  = c(0.25, 0.5, 0.75, 1)),
    # relative increment of clinical immunity following relapse w.r.t. reinfection
    tibble(a_clin_rel = c(0.25, 0.5, 0.75, 1))
    
  )
  
}
