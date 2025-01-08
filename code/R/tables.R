# tab_eligibility_criteria <- function () {
#   
#   tab_name <- "eligibility_criteria"
#   
#   tibble::tribble(
#     # ----------------------------------------------------------------------------------------------------------------------------- #
#     ~`Drug regimen`,               ~`Eligibility criteria`,
#     # ----------------------------------------------------------------------------------------------------------------------------- #
#     "No 8-AQ",                     "none",
#     "PQ (3.5 mg/kg over 14 days)", "aged > 6m, G6PD activity > 30%, not pregnant, not breastfeeding",
#     "PQ (3.5 mg/kg over 7 days)",  "aged > 6m, G6PD activity > 30%, not pregnant, not breastfeeding",
#     "PQ (7 mg/kg over 14 days)",   "aged > 6m, G6PD activity > 30%, not pregnant, not breastfeeding",
#     "PQ (7 mg/kg over 7 days)",    "aged > 6m, G6PD activity > 70%, not pregnant, not breastfeeding",
#     "TQ (5 mg/kg single dose)",    "aged > 2y *, G6PD activity > 70%, not pregnant, not breastfeeding",
#     "TQ (7.5 mg/kg single dose)",  "aged > 2y *, G6PD activity > 70%, not pregnant, not breastfeeding",
#     "Magic bullet (single dose)",  "none (assumption)"
#   ) %>%
#     # print to file
#     write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
#   
# }


tab_drug_regimen_overview <- function () {
  
  tab_name <- "drug_regimen_overview"
  
  tibble::tribble(
    # ----------------------------------------------------------------------------------------------------------------------------- #
    ~`8-AQ regimen`,               ~`8-AQ`,        ~nd,  ~dd,    ~td,   ~td60, ~dec,                                        ~gec,
    # ----------------------------------------------------------------------------------------------------------------------------- #
    "No 8-AQ",                     "none",         "0",  "N/A",  "N/A", "N/A", "none",                                            "none",
    "PQ (3.5 mg/kg over 14 days)", "PQ",           "14", "0.25", "3.5", "210", "aged >6 months, not pregnant, not breastfeeding", "G6PD activity >30%",
    "PQ (3.5 mg/kg over 7 days)",  "PQ",           "7",  "0.5",  "3.5", "210", "aged >6 months, not pregnant, not breastfeeding", "G6PD activity >30%",
    "PQ (7 mg/kg over 14 days)",   "PQ",           "14", "0.5",  "7",   "420", "aged >6 months, not pregnant, not breastfeeding", "G6PD activity >30%",
    "PQ (7 mg/kg over 7 days)",    "PQ",           "7",  "1",    "7",   "420", "aged >6 months, not pregnant, not breastfeeding", "G6PD activity >70%",
    "TQ (5 mg/kg single dose)",    "TQ",           "1",  "5",    "5",   "300", "aged >2 years,  not pregnant, not breastfeeding", "G6PD activity >70%",
    "TQ (7.5 mg/kg single dose)",  "TQ",           "1",  "7.5",  "7.5", "450", "aged >2 years,  not pregnant, not breastfeeding", "G6PD activity >70%",
    "Magic bullet (single dose)",  "hypothetical", "1",  "N/A",  "N/A", "N/A", "none",                                            "none"
  ) %>%
    rename(
      `Number of days`                      = nd,
      `Daily dose (mg/kg)`                  = dd,
      `Total dose (mg/kg)`                  = td,
      `Total dose for a 60-kg patient (mg)` = td60,
      `Demographic eligibility criteria`    = dec,
      `G6PD eligibility criteria*`          = gec
    ) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_estimates_rrela_rbite <- function (taskname) {
  
  tab_name <- paste("ABC-estimates-rate", taskname, sep = "-")
  
  
  # transform
  dt <-
    lapply(c("A", "B", "C"), function (step) {
      
      read_estimates(step, "MCMC", taskname) %>%
        compute_proportion_of_relapses() %>%
        pivot_and_rescale() %>%
        compute_summary_stats() %>%
        mutate(source = step)
      
    }) %>%
    bind_rows()
  
  
  # subset other fitted parameters
  dt %>%
    filter(grepl("^rbite|^rrela|^prela", name)) %>%
    select(name3, line1, line2, line3, source) %>%
    pivot_wider(names_from = line1, values_from = line3) %>%
    rename(Location = line2) %>%
    # clean
    mutate(
      `Fitted to` = factor(source, names(label_source), label_source),
      source      = NULL
    ) %>%
    arrange(name3) %>%
    select(-name3) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_proportion_eligible_for_8AQ <- function () {
  
  tab_name <- "proportion_eligible_for_8AQ"
  
  dt_default %>%
    transmute(
      
      # demography
      demog__less6m = prop_aged_less(PQ__lowage, age_mean, age_max),
      demog__less2y = prop_aged_less(TQ__lowage, age_mean, age_max),
      demog__fertile_age = prop_aged_less(age_preg_max, age_mean, age_max) - prop_aged_less(age_preg_min, age_mean, age_max),
      
      
      # females
      f__pregnant_or_breastfeeding = P_preg * demog__fertile_age,
      f__over6m_not_pregnant       = 1 - demog__less6m - f__pregnant_or_breastfeeding,
      f__over2y_not_pregnant       = 1 - demog__less2y - f__pregnant_or_breastfeeding,
      
      f__G6PD_over30 = 1 - prop_G6PD_activity_less(sex = "female", threshold = 3, G6PD_prev, mu_G6PD_nor, sig_G6PD_nor, mu_G6PD_het, sig_G6PD_het, mu_G6PD_def, sig_G6PD_def),
      f__G6PD_over70 = 1 - prop_G6PD_activity_less(sex = "female", threshold = 7, G6PD_prev, mu_G6PD_nor, sig_G6PD_nor, mu_G6PD_het, sig_G6PD_het, mu_G6PD_def, sig_G6PD_def),
      
      
      # males
      m__over6m = 1 - demog__less6m,
      m__over2y = 1 - demog__less2y,
      
      m__G6PD_over30 = 1 - prop_G6PD_activity_less(sex = "male", threshold = 3, G6PD_prev, mu_G6PD_nor, sig_G6PD_nor, mu_G6PD_het, sig_G6PD_het, mu_G6PD_def, sig_G6PD_def),
      m__G6PD_over70 = 1 - prop_G6PD_activity_less(sex = "male", threshold = 7, G6PD_prev, mu_G6PD_nor, sig_G6PD_nor, mu_G6PD_het, sig_G6PD_het, mu_G6PD_def, sig_G6PD_def),
      
      
      # eligible proportions
      elig_PQ      = prop_female * (f__over6m_not_pregnant * f__G6PD_over30) + (1 - prop_female) * (m__over6m * m__G6PD_over30),
      elig_PQrestr = prop_female * (f__over6m_not_pregnant * f__G6PD_over70) + (1 - prop_female) * (m__over6m * m__G6PD_over70),
      elig_TQ      = prop_female * (f__over2y_not_pregnant * f__G6PD_over70) + (1 - prop_female) * (m__over2y * m__G6PD_over70),
      
    ) %>%
    # writ to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


# add hypnozoiticidal efficacy
add_hypeff <- function (row_) {
  
  taskname <- "AON"
  
  dt_hypeff <-
    lapply(c("A", "B", "C"), function (step) {
      
      read_estimates(step, "MCMC", taskname) %>%
        pivot_and_rescale() %>%
        compute_summary_stats() %>%
        # clean
        filter(grepl("^eff__", name)) %>%
        select(name3, m, q1, q2, line3)
      
    }) %>%
    bind_rows() %>%
    rename(name = name3, CI = line3)
  
  
  m   <- dt_hypeff %>% select(name, m)   %>% pivot_wider(names_from = name, values_from = m)
  q1  <- dt_hypeff %>% select(name, q1)  %>% pivot_wider(names_from = name, values_from = q1)
  q2  <- dt_hypeff %>% select(name, q2)  %>% pivot_wider(names_from = name, values_from = q2)
  CI  <- dt_hypeff %>% select(name, CI)  %>% pivot_wider(names_from = name, values_from = CI)
  
  case_when(
    grepl("No_8", row_) ~ tibble(he_m = 0,               he_q1 = 0,                he_q2 = 0,                he = paste("0 [0, 0] ^a^")),
    grepl("AQ_i", row_) ~ tibble(he_m = 100,             he_q1 = 100,              he_q2 = 100,              he = paste("100 [100, 100] ^a^")),
    grepl("PQ_l", row_) ~ tibble(he_m = m$PQ_lowdose,    he_q1 = q1$PQ_lowdose,    he_q2 = q2$PQ_lowdose,    he = paste(CI$PQ_lowdose,    "^b^")),
    grepl("PQ_h", row_) ~ tibble(he_m = m$PQ_highdose14, he_q1 = q1$PQ_highdose14, he_q2 = q2$PQ_highdose14, he = paste(CI$PQ_highdose14, "^c^")),
    grepl("PQr",  row_) ~ tibble(he_m = m$PQ_highdose7,  he_q1 = q1$PQ_highdose7,  he_q2 = q2$PQ_highdose7,  he = paste(CI$PQ_highdose7,  "^c^")),
    grepl("TQ_l", row_) ~ tibble(he_m = m$TQ_lowdose,    he_q1 = q1$TQ_lowdose,    he_q2 = q2$TQ_lowdose,    he = paste(CI$TQ_lowdose,    "^d^")),
    grepl("TQ_h", row_) ~ tibble(he_m = m$TQ_highdose,   he_q1 = q1$TQ_highdose,   he_q2 = q2$TQ_highdose,   he = paste(CI$TQ_highdose,   "^d^")),
    grepl("d1",   row_) ~ tibble(he = paste0("Patients who")),
    grepl("d2",   row_) ~ tibble(he = paste0("- are eligible for the 8-AQ regimen")),
    grepl("d3",   row_) ~ tibble(he = paste0("- perfectly adhere")),
    grepl("d4",   row_) ~ tibble(he = paste0("- are CYP2D6-normal")),
    grepl("n1",   row_) ~ tibble(he = paste0("Inference framework:")),
    grepl("n2",   row_) ~ tibble(he = paste0("^a^ Assumption")),
    grepl("n3",   row_) ~ tibble(he = paste0("^b^ Fitted to meta-analysis of patient-level trial data [@Commons2023]")),
    grepl("n4",   row_) ~ tibble(he = paste0("^c^ Fitted to patient-level clinical trial data [@Taylor2019a]")),
    grepl("n5",   row_) ~ tibble(he = paste0("^d^ Fitted to meta-analysis of patient-level trial data [@Watson2022]"))
  )
  
}

# add median effectiveness under trial conditions
add_eff_tc <- function (row_, he_m, he_q1, he_q2) {
  
  case_when(
    grepl("^PQ", row_) ~ tibble(
      etc_m  = he_m  * (1 - dt_default$CYP2D6_prev),
      etc_q1 = he_q1 * (1 - dt_default$CYP2D6_prev),
      etc_q2 = he_q2 * (1 - dt_default$CYP2D6_prev),
      etc    = paste0(round(etc_m, 1), " [", round(etc_q1, 1), ", ", round(etc_q2, 1), "] ^e^")
    ),
    grepl("_",   row_) ~ tibble(
      etc_m  = he_m,
      etc_q1 = he_q1,
      etc_q2 = he_q2,
      etc    = paste0(round(etc_m, 1), " [", round(etc_q1, 1), ", ", round(etc_q2, 1), "]")
    ),
    grepl("d1",  row_) ~ tibble(etc = "Patients who"),
    grepl("d2",  row_) ~ tibble(etc = "- are eligible for the 8-AQ regimen"),
    grepl("d3",  row_) ~ tibble(etc = "- perfectly adhere"),
    grepl("d4",  row_) ~ tibble(etc = NA_character_),
    grepl("n1",  row_) ~ tibble(etc = "^e^ patients who are low CYP2D6 metabolisers are unable to process PQ"),
    grepl("n2",  row_) ~ tibble(etc = NA_character_),
    grepl("n3",  row_) ~ tibble(etc = NA_character_),
    grepl("n4",  row_) ~ tibble(etc = NA_character_),
    grepl("n5",  row_) ~ tibble(etc = NA_character_),
  )
  
}

# add median effectiveness under operational conditions
add_eff_oc <- function (dt_elig, row_, etc_m, etc_q1, etc_q2) {
  
  
  case_when(
    grepl("PQr",   row_) ~ tibble(
      eoc_m  = etc_m  * dt_elig$elig_PQrestr * dt_default$PQ_7__adhere,
      eoc_q1 = etc_q1 * dt_elig$elig_PQrestr * dt_default$PQ_7__adhere,
      eoc_q2 = etc_q2 * dt_elig$elig_PQrestr * dt_default$PQ_7__adhere,
      eoc    = paste0(round(eoc_m, 1), " [", round(eoc_q1, 1), ", ", round(eoc_q2, 1), "]")
    ),
    grepl("7$",  row_) ~ tibble(
      eoc_m  = etc_m  * dt_elig$elig_PQ * dt_default$PQ_7__adhere,
      eoc_q1 = etc_q1 * dt_elig$elig_PQ * dt_default$PQ_7__adhere,
      eoc_q2 = etc_q2 * dt_elig$elig_PQ * dt_default$PQ_7__adhere,
      eoc    = paste0(round(eoc_m, 1), " [", round(eoc_q1, 1), ", ", round(eoc_q2, 1), "]")
    ),
    grepl("14$", row_) ~ tibble(
      eoc_m  = etc_m  * dt_elig$elig_PQ * dt_default$PQ_14__adhere,
      eoc_q1 = etc_q1 * dt_elig$elig_PQ * dt_default$PQ_14__adhere,
      eoc_q2 = etc_q2 * dt_elig$elig_PQ * dt_default$PQ_14__adhere,
      eoc    = paste0(round(eoc_m, 1), " [", round(eoc_q1, 1), ", ", round(eoc_q2, 1), "]")
    ),
    grepl("TQ",   row_) ~ tibble(
      eoc_m  = etc_m  * dt_elig$elig_TQ,
      eoc_q1 = etc_q1 * dt_elig$elig_TQ,
      eoc_q2 = etc_q2 * dt_elig$elig_TQ,
      eoc    = paste0(round(eoc_m, 1), " [", round(eoc_q1, 1), ", ", round(eoc_q2, 1), "]")
    ),
    grepl("_",   row_) ~ tibble(
      eoc_m  = etc_m,
      eoc_q1 = etc_q1,
      eoc_q2 = etc_q2,
      eoc    = paste0(round(eoc_m, 1), " [", round(eoc_q1, 1), ", ", round(eoc_q2, 1), "]")
    ),
    grepl("d1",  row_) ~ tibble(eoc = "Patients who"),
    grepl("d2",  row_) ~ tibble(eoc = "- seek treatment"),
    grepl("d3",  row_) ~ tibble(eoc = NA_character_),
    grepl("d4",  row_) ~ tibble(eoc = NA_character_),
    grepl("n1",  row_) ~ tibble(eoc = NA_character_),
    grepl("n2",  row_) ~ tibble(eoc = NA_character_),
    grepl("n3",  row_) ~ tibble(eoc = NA_character_),
    grepl("n4",  row_) ~ tibble(eoc = NA_character_),
    grepl("n5",  row_) ~ tibble(eoc = NA_character_)
  )
  
}

# add median effectiveness under operational conditions adding secondary treatment
add_eff_oc2 <- function (dt_elig, row_, eoc_m, eoc_q1, eoc_q2, eoc, v) {
  
  case_when(
    grepl("PQr", row_) ~ tibble(
      eoc2_m  = eoc_m  + ( (dt_elig$elig_PQ - dt_elig$elig_PQrestr) * v$etc_m  * dt_default$PQ_7__adhere ),
      eoc2_q1 = eoc_q1 + ( (dt_elig$elig_PQ - dt_elig$elig_PQrestr) * v$etc_q1 * dt_default$PQ_7__adhere ),
      eoc2_q2 = eoc_q2 + ( (dt_elig$elig_PQ - dt_elig$elig_PQrestr) * v$etc_q2 * dt_default$PQ_7__adhere ),
      eoc2    = paste0(round(eoc2_m, 1), " [", round(eoc2_q1, 1), ", ", round(eoc2_q2, 1), "] ^f^")
    ),
    grepl("TQ", row_) ~ tibble(
      eoc2_m  = eoc_m  + ( (dt_elig$elig_PQ - dt_elig$elig_TQ) * v$etc_m  * dt_default$PQ_7__adhere ),
      eoc2_q1 = eoc_q1 + ( (dt_elig$elig_PQ - dt_elig$elig_TQ) * v$etc_q1 * dt_default$PQ_7__adhere ),
      eoc2_q2 = eoc_q2 + ( (dt_elig$elig_PQ - dt_elig$elig_TQ) * v$etc_q2 * dt_default$PQ_7__adhere ),
      eoc2    = paste0(round(eoc2_m, 1), " [", round(eoc2_q1, 1), ", ", round(eoc2_q2, 1), "] ^f^")
    ),
    grepl("n1",  row_) ~ tibble(eoc2 = "Alternative 8-AQ regimen:"),
    grepl("n2",  row_) ~ tibble(
      eoc2 = paste0("^f^ patients who are not eligible for the use of TQ or PQ (7 mg/kg over 7 days), ",
                    "are prescribed PQ (3.5 mg/kg over 7 days) provided they have a G6PD activity >30%, ",
                    "are not pregnant or breastfeeding, and are 6 months or older")),
    grepl("n3",  row_) ~ tibble(eoc2 = NA_character_),
    grepl("n4",  row_) ~ tibble(eoc2 = NA_character_),
    grepl("n5",  row_) ~ tibble(eoc2 = NA_character_),
    .default = tibble(
      eoc2_m  = eoc_m,
      eoc2_q1 = eoc_q1,
      eoc2_q2 = eoc_q2,
      eoc2    = eoc)
  )
  
}


tab_estimates_efficacy_effectiveness <- function (taskname) {
  
  tab_name <- paste("ABC-estimates-eff", taskname, sep = "-")
  
  
  # read in proportion of IBM population eligible for 8-AQ regimens
  dt_elig <- read_csv(file.path(dir_tab, "proportion_eligible_for_8AQ.csv"))
  
  
  # initialise table
  dt <-
    bind_rows(
      # Upper group containing estimates
      tibble(
        group = "8-AQ regimen",
        row   = label_drug_regimen %>% filter(!is.na(colour)) %>% pull(label),
        row_  = label_drug_regimen %>% filter(!is.na(colour)) %>% pull(name)
      ),
      # Bottom group containing descriptions
      tibble(
        group = "Description",
        row   = c("Denominator", rep("", 3), "Notes", rep("", 4)),
        row_  = c(paste0("d", 1:4),          paste0("n", 1:5))
      )
    ) %>%
    # populate table
    mutate(
      # add hypnozoiticidal efficacy estimates
      add_hypeff(row_),
      # add median effectiveness under trial conditions
      add_eff_tc(row_, he_m, he_q1, he_q2),
      # add median effectiveness under operational conditions
      add_eff_oc(dt_elig, row_, etc_m, etc_q1, etc_q2)
    )
  
  
  # get key estimates for secondary treatment PQ_lowdose7
  v_PQld7 <- dt %>% filter(row_ == "PQ_lowdose7") %>% select(etc_m, etc_q1, etc_q2)
  
  
  dt %>%
    # add median effectiveness under operational conditions adding secondary treatment
    mutate(add_eff_oc2(dt_elig, row_, eoc_m, eoc_q1, eoc_q2, eoc, v_PQld7)) %>%
    # writ to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_trial_data_overview <- function () {
  
  tab_name <- "ABC-trial_data_overview"
  
  dtA <-
    read_csv(file.path(dirX_inp("A"), "data.csv")) %>%
    group_by(Drug_regimen) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Drug_regimen, values_from = n)
  
  dtB <-
    read_csv(file.path(dirX_inp("B"), "data.csv")) %>%
    pivot_wider(names_from = Drug_regimen, values_from = nr_enrolled, id_cols = !everything())
  
  dtC <-
    read_csv(file.path(dirX_inp("C"), "data.csv")) %>%
    pivot_wider(names_from = Drug_regimen, values_from = nr_enrolled, id_cols = !everything())
  
  tibble::tribble(
    # ----------------------------------------------------------------------------------------------------------------------------- #
    ~empty,                     ~`IMPROV trial [@Taylor2019a]`,              ~`Meta-analysis [@Commons2023]`,            ~`Meta-analysis [@Watson2022]`,
    # ----------------------------------------------------------------------------------------------------------------------------- #
    "PQ supervised",            "Fully (100%)",                              "Fully (63.4%), partially (33.5%)",         "Fully (32.2%), partially (67.8%)",
    "Partner blood-stage drug", "CQ, except DP in Indonesia",                "CQ (67.1%), DP (23.1%), other ACT (9.8%)", "CQ",
    "Follow-up period",         "365 days",                                  "Varied, minimum 42 days",                  "180 days",
    "Recurrence endpoint",      "First symptomatic *P. vivax* recurrence",   "First *P. vivax* recurrence",              "First *P. vivax* recurrence",
    "Locations",                "Afghanistan, Ethiopia, Indonesia, Vietnam", "Asia-Pacific, Americas, Africa",           "Asia-Pacific, Americas, Africa",
    "Trial arms",
    paste0(
      "Placebo (n=", dtA$No_8AQs, "), ",
      "PQ 7 mg/kg over 7 days (n=", dtA$PQ_highdose7, "), ",
      "PQ 7 mg/kg over 14 days (n=", dtA$PQ_highdose14, ")"
    ),
    paste0(
      "Control* (n=", dtB$No_8AQs, "), ",
      "PQ 3.5 mg/kg (n=", dtB$PQ_lowdose, "), ",
      "PQ 7 mg/kg (n=", dtB$PQ_highdose, ")"
    ),
    paste0(
      "Placebo (n=", dtC$No_8AQs, "), ",
      "PQ 3.5 mg/kg over 14 days (n=", dtC$PQ_lowdose14, "), ",
      "TQ 5 mg/kg single dose (n=", dtC$TQ_lowdose, "), ",
      "TQ 7.5 mg/kg single dose (n=", dtC$TQ_highdose, ")"
    )
  ) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_KM_helper <- function (dt, limit) {
  
  colname <- paste("Day", limit)
  
  dt %>%
    filter(Day < limit) %>%
    group_by(Location, Drug_regimen) %>%
    slice_max(Day) %>%
    transmute(nr_survived) %>%
    ungroup() %>%
    rename(!!colname := nr_survived)
  
}


tab_KM <- function () {
  
  tab_name <- "A-KM"
  
  dt <-
    read_data_A(step = "A") %>%
    group_by(Location, Drug_regimen) %>%
    mutate(
      temp = max(nr_enrolled, na.rm = TRUE),
      nr_survived = if_else(is.na(nr_survived), temp, nr_survived)
    ) %>%
    ungroup()
  
  plyr::join_all(
    list(
      tab_KM_helper(dt, limit = 42),
      tab_KM_helper(dt, limit = 86),
      tab_KM_helper(dt, limit = 176),
      tab_KM_helper(dt, limit = 267),
      tab_KM_helper(dt, limit = 356)
    ),
    by = c("Location", "Drug_regimen")
  ) %>%
    mutate(
      Location     = factor(Location,     label_location$name,     label_location$label2),
      Drug_regimen = factor(Drug_regimen, label_drug_regimen$name, label_drug_regimen$label)
    ) %>%
    rename(`8-AQ regimen` = Drug_regimen) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_trial_data_A <- function () {
  
  tab_name <- "A-trial_data"
  
  read_csv(file.path(dirX_inp("A"), "data.csv")) %>%
    group_by(Location, Drug_regimen) %>%
    summarise(
      `Nr enrolled`            = n(),
      `Nr recurred by day 365` = sum(fail)
    ) %>%
    ungroup() %>%
    mutate(
      Location     = factor(Location,     label_location$name,     label_location$label2),
      Drug_regimen = factor(Drug_regimen, label_drug_regimen$name, label_drug_regimen$label)
    ) %>%
    rename(`8-AQ regimen` = Drug_regimen) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_trial_data_B <- function () {
  
  tab_name <- "B-trial_data"
  
  dt_B <- read_csv(file.path(dirX_inp("B"), "data.csv"))
  
  dt_B %>%
    mutate(
      m  = 100 * m,
      q1 = 100 * q1,
      q2 = 100 * q2,
      Drug_regimen = factor(Drug_regimen, label_drug_regimen$name, label_drug_regimen$label),
      prop_recurrence = paste0(m, " [", q1, ", ", q2, "]")
    ) %>%
    rename(
      `8-AQ regimen` = Drug_regimen,
      `Nr enrolled`  = nr_enrolled
    ) %>%
    rename_with(
      ~ paste0("Proportion recurred by day ", unique(dt_B$Day), " (%)"),
      prop_recurrence
    ) %>%
    select(!c(m, q1, q2, Day)) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_trial_data_C <- function () {
  
  tab_name <- "C-trial_data"
  
  dt_C <- read_csv(file.path(dirX_inp("C"), "data.csv"))
  
  dt_C %>%
    mutate(
      Drug_regimen = factor(Drug_regimen, label_drug_regimen$name, label_drug_regimen$label)
    ) %>%
    rename(
      `8-AQ regimen` = Drug_regimen,
      `Nr enrolled`  = nr_enrolled
    ) %>%
    rename_with(
      ~ paste0("Nr recurred by day ", unique(dt_C$Day)),
      nr_recurrences
    ) %>%
    select(!Day) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_varying_recurrence_model <- function () {
  
  tab_name <- "ABC-varying"
  
  dtA <- read_csv(file.path(dirX_inp("A"), "varying.csv")) %>% select(!c(AON, country))
  dtC <- read_csv(file.path(dirX_inp("C"), "varying.csv")) %>% select(!AON)
  
  # get default values
  default <- 
    full_join(
      dtA %>% filter(taskname == "AON"),
      dtC %>% filter(taskname == "AON")
    ) %>%
    select(!taskname)
  
  # get alternative values
  varying <- bind_rows(dtA, dtC) %>% select(!taskname) %>% unique()
  
  varying <-
    lapply(
      # loop over column names
      names(varying),
      function (x) {
        
        varying[x] %>%
          # get column of unique values
          unique() %>%
          # arrange ascendingly
          arrange(across(everything())) %>%
          # filter values that are NOT NA and NOT default values
          filter(.data[[x]] != default[[x]]) %>%
          # make into a string
          summarise(across(everything(), ~ paste(.x, collapse = ", ")))
        
      }) %>%
    bind_cols()
  
  
  # Parameter description, notation and reference for default value
  dt_descr <-
    tribble(
      ~name,         ~Description,                                     ~Parameter,    ~`Reference for default value`,
      "CYP2D6_prev", "Prevalence of CYP2D6 low metabolisers",          "$c$",         "@Koopmans2021",
      "d_gamma",     "Mean duration of hypnozoite carriage (in days)", "$1/\\gamma$", "@White2014a",
      "proph_CQ",    "Duration of CQ prophylaxis (in days)",           "$d$",         "@Karunajeewa2008",
      "proph_DP",    "Duration of DP prophylaxis (in days)",           "$d$",         "@Karunajeewa2008",
      "proph_TQl",   "Duration of TQ 5 mg/kg prophylaxis (in days)",   "$d$",         "@Lacerda2019",
      "proph_TQh",   "Duration of TQ 7.5 mg/kg prophylaxis (in days)", "$d$",         "Assumption"
    )
  
  
  # combine
  inner_join(
    default %>% pivot_longer(cols = everything(), values_to = "Default value"),
    varying %>% pivot_longer(cols = everything(), values_to = "Alternative values"),
    by = join_by(name)
  ) %>%
    inner_join(
      dt_descr,
      by = join_by(name)
    ) %>%
    # clean
    select(Description, Parameter, `Default value`, `Alternative values`, `Reference for default value`) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_varying_IBM <- function () {
  
  tab_name <- "D-varying"
  
  default <-
    default_parameters(print_to_file = T) %>%
    select(!References) %>%
    rename(`Default value` = Value)
  
  varying <-
    lapply(
      # loop over column names
      names(dt_varying),
      function (x) {
        
        dt_varying[x] %>%
          # get column of unique values
          unique() %>%
          # arrange ascendingly
          arrange(across(everything())) %>%
          # filter values that are NOT NA and NOT default values
          filter(.data[[x]] != dt_default[[x]]) %>%
          # make into a string
          summarise(across(everything(), ~ paste(.x, collapse = ", ")))
        
      }) %>%
    bind_cols() %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Alternative values"
    )
  
  right_join(default, varying) %>%
    select(!Parameter) %>%
    # remove occupational exposure parameters
    filter(!grepl("occup", Description)) %>%
    relocate(Description) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_posterior_draws_IBM <- function () {
  
  tab_name <- "D-posterior_draws"
  
  default <-
    default_parameters(print_to_file = T) %>%
    select(Parameter, Description)
  
  posterior <-
    # read in the posterior parameter draws obtained from the calibration of the 
    # PvIBM by [White et al](http://dx.doi.org/10.1038/s41467-018-05860-8)
    read_csv(file.path(dir_raw, "PvMod_posterior_draws.csv")) %>%
    # rename
    rename(
      A_PCR_50pc = A_d_PCR_50pc,
      K_PCR      = K_d_PCR,
      P_mat      = P_MI,
      d_mat      = d_MI
    ) %>%
    # summarise
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(
      m   = median(value),
      q1  = quantile(value, probs = 0.025),
      q2  = quantile(value, probs = 0.975)
    ) %>%
    ungroup() %>%
    # CrI
    mutate(
      nd = if_else(q1 > 9, 2, 3),
      CrI = paste0(round(m, nd), " [", round(q1, nd), ", ", round(q2, nd), "]")
    ) %>%
    # clean
    select(name, CrI) %>%
    pivot_wider(names_from = name, values_from = CrI) %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Median and 95% CI")
  
  right_join(default, posterior) %>%
    select(!Parameter) %>%
    relocate(Description) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}


tab_fixed_IBM <- function () {
  
  tab_name <- "D-fixed"
  
  default_parameters(print_to_file = T) %>%
    filter(!is.na(Value) & References != "") %>%
    select(Description, Value, References) %>%
    mutate(Value = Value %>% signif(3) %>% as.character()) %>%
    # print to file
    write_csv(file = file.path(dir_tab, paste0(tab_name, ".csv")))
  
}
