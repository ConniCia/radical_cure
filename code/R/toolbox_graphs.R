
# colours ----------------------------------------------------------------------

# colours               # Drug regimens     # Treatment courses
my.black  <- "black"    # No_8AQs           # sum_new_CQ
my.red    <- "#e41a1c"  # PQ_lowdose14      # sum_G6PD_tests
my.yellow <- "#ffd92f"  # PQ_lowdose7       #
my.brown  <- "#a65628"  # PQ_highdose14     #
my.blue   <- "#377eb8"  # PQrestr_highdose7 # sum_new_PQ
my.pink   <- "#f781bf"  # PQ_highdose7      #
my.purple <- "#984ea3"  # TQ_lowdose        #
my.orange <- "#ff7f00"  # TQ_highdose       # sum_new_TQ
my.green  <- "#4daf4a"  # AQ_ideal          # sum_new_AQ


# labels -----------------------------------------------------------------------

label_drug_regimen <-
  tribble(
    # -------------------------------------------------------------------------------------- #
    ~name,               ~colour,   ~label,                        ~label2,
    # -------------------------------------------------------------------------------------- #
    "No_8AQs",           my.black,  "No 8-AQs",                    "No 8-AQs",
    "PQ_lowdose",        NA,        "PQ (3.5 mg/kg)",              "PQ (3.5 mg/kg)",
    "PQ_lowdose14",      my.red,    "PQ (3.5 mg/kg over 14 days)", "PQ (3.5 mg/kg\nover 14 days)",
    "PQ_lowdose7",       my.yellow, "PQ (3.5 mg/kg over 7 days)",  "PQ (3.5 mg/kg\nover 7 days)",
    "PQ_highdose",       NA,        "PQ (7 mg/kg)",                "PQ (7 mg/kg)",
    "PQ_highdose14",     my.brown,  "PQ (7 mg/kg over 14 days)",   "PQ (7 mg/kg\nover 14 days)",
    "PQ_highdose7",      NA,        "PQ (7 mg/kg over 7 days)",    "PQ (7 mg/kg\nover 7 days)",
    "PQrestr_highdose7", my.blue,   "PQ (7 mg/kg over 7 days)",    "PQ (7 mg/kg\nover 7 days)",
    "TQ_lowdose",        my.purple, "TQ (5 mg/kg single dose)",    "TQ (5 mg/kg\nsingle dose)",
    "TQ_highdose",       my.orange, "TQ (7.5 mg/kg single dose)",  "TQ (7.5 mg/kg\nsingle dose)",
    "AQ_ideal",          my.green,  "Magic bullet (single dose)",  "Magic bullet\n(single dose)"
  )

label_treatment_courses <-
  tribble(
    # -------------------------------------- #
    ~name,            ~colour,   ~label,
    # -------------------------------------- #
    "sum_new_CQ",     my.black,  "CQ courses",
    "sum_G6PD_tests", my.red,    "G6PD tests",
    "sum_new_PQ",     my.blue,   "PQ courses",
    "sum_new_TQ",     my.orange, "TQ courses",
    "sum_new_AQ",     my.green,  "Magic bullet courses"
  )

label_relapses_reinfections <-
  tribble(
    # -------------------------------------- #
    ~name,                ~colour,  ~label,
    # -------------------------------------- #
    "sum_new_CQ_bite",    my.black, "Reinfection",
    "sum_new_CQ_relapse", my.red,   "Relapse"
  )

label_fitted <- c(
  eff   = "Hypnozoiticidal efficacy (%)",
  rbite = "Yearly reinfection rate",
  rrela = "Yearly relapse rate",
  prela = "Proportion of recurrences due to relapses (%)",
  ll    = "Log-likelihood"
)

label_source <- c(
  A = "Patient-level trial data [@Taylor2019a]",
  B = "Meta-analysis of patient-level trial data [@Commons2023]",
  C = "Meta-analysis of patient-level trial data [@Watson2022]"
)

label_location <- 
  tribble(
    # ------------------------------------------------------------------------ #
    ~name,   ~label,                             ~label2,
    # ------------------------------------------------------------------------ #
    "AF001", "Jalalabad,\nAfghanistan",          "Jalalabad, Afghanistan",
    "AF008", "Laghman,\nAfghanistan",            "Laghman, Afghanistan",
    "ET001", "Arba Minch,\nEthiopia",            "Arba Minch, Ethiopia",
    "ET002", "Metahara,\nEthiopia",              "Metahara, Ethiopia",
    "ID004", "Hanura,\nIndonesia",               "Hanura, Indonesia",
    "ID005", "Tanjung Leidong,\nIndonesia",      "Tanjung Leidong, Indonesia",
    "VN001", "Dak O and\nBu Gia Map,\nVietnam",  "Dak O and Bu Gia Map, Vietnam",
    "VN002", "Krong Pa,\nVietnam",               "Krong Pa, Vietnam",
    "B",     "Various locations",                "Various locations",
    "C",     "Various locations ",               "Various locations "
  )

label_conditions <- c(
  tc = "Trial conditions",
  oc = "Operational conditions"
)

label_conditions2 <- c(
  tc = "Trial\nconditions",
  oc = "Operational\nconditions"
)

label_pars_varying <- c(
  a_clin_rel  = "Increment of clinical immunity\nfollowing a relapse",
  a_par_rel   = "Increment of anti-parasite\nimmunity following a relapse",
  CM_cover    = "Case management\ncoverage",
  CYP2D6_prev = "CYP2D6 prevalence",
  d_ff        = "Relapse periodicity (days)",
  G6PD_prev   = "G6PD prevalence",
  phi_D_rel   = "Relative proportion of relapses\nthat cause symptomatic infection",
  phi_LM_rel  = "Relative proportion of relapses\nthat cause microscopic infection",
  PQ__adhere  = "PQ adherence",
  seasonality = "Seasonality",
  yEIR        = "Yearly EIR at equilibrium"
)

label_name <- c(
  effsize       = "Reduction in PCR-prevalence\nafter 5 years (%)",
  cases_averted = "Averted clinical cases over\n5 years (per 1 000 population)",
  mean_prev_PCR = "PCR-prevalence\npre intervention (%)"
)

order_fitted_pars <- c(
  paste0("rrela__", label_location$name),
  paste0("rbite__", label_location$name),
  paste0("prela__", label_location$name),
  paste0("eff__",   label_drug_regimen$name),
  "ll"
)

label_transmission_intensity <- c(
  "Low transmission intensity\n(~2% baseline PCR-prevalence)",
  "Moderate transmission intensity\n(~10% baseline PCR-prevalence)",
  "High transmission intensity\n(~35% baseline PCR-prevalence)"
)

label_transmission_intensity2 <- c(
  "Low transmission intensity\n(~2% PCR-prevalence)",
  "Moderate transmission intensity\n(~10% PCR-prevalence)",
  "High transmission intensity\n(~35% PCR-prevalence)"
)


global_labeller <- 
  labeller(
    conditions   = label_conditions,
    Conditions   = label_conditions2,
    drug_regimen = tibble::deframe(label_drug_regimen %>% select(name, label)),
    Drug_regimen = tibble::deframe(label_drug_regimen %>% select(name, label2)),
    Location     = tibble::deframe(label_location %>% select(name, label)),
    name         = label_name,
    pars_varying = label_pars_varying,
    variable = c(
      smooth_prev_PCR = "Median PCR-prevalence (%)",
      smooth_new_T    = "Median clinical incidence\n(per 1 000 population)"
    ),
    yEIR         = function (str) {paste(str, "infectious bites per person per year at equilibrium")},
    yEIR2        = function (str) {paste(str, "infectious bites per person per year")},
    yeir         = function (str) {x <- label_transmission_intensity;  names(x) <- str; x},
    yeir2        = function (str) {x <- label_transmission_intensity2; names(x) <- str; x}
  )


# discrete colour/fill scales --------------------------------------------------

my.scale_drug_regimen <- function (aesthetics, ...) {
  scale_colour_manual(
    aesthetics = aesthetics,
    name       = "8-AQ regimen",
    values     = tibble::deframe(label_drug_regimen %>% select(name, colour)),
    labels     = tibble::deframe(label_drug_regimen %>% select(name, label)),
    ...
  )
}

my.scale_treatment_courses <- function (aesthetics, ...) {
  scale_colour_manual(
    aesthetics = aesthetics,
    name       = "Diagnostic tests and treatment courses",
    values     = tibble::deframe(label_treatment_courses %>% select(name, colour)),
    labels     = tibble::deframe(label_treatment_courses %>% select(name, label)),
    ...
  )
}

my.scale_relapses_reinfections <- function (aesthetics, ...) {
  scale_colour_manual(
    aesthetics = aesthetics,
    name       = "Number of recurrences",
    values     = tibble::deframe(label_relapses_reinfections %>% select(name, colour)),
    labels     = tibble::deframe(label_relapses_reinfections %>% select(name, label)),
    ...
  )
}
