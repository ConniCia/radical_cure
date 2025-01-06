clean_D <- function (current_dir) {
  
  # Read in values of looped parameters corresponding to current directory
  row_varying <-
    read_delim(file.path(current_dir, "row_varying.tsv"), delim = " ") %>%
    select(-c("sample_nr", starts_with("path_")))
  
  
  # Define column names
  cols <-
    c(
      # 1: time in days
      "day",
      # 2-7: number of domestic mosquitoes in each compartment (per individual????)
      "EL_M_dom", "LL_M_dom", "P_M_dom", "S_M_dom", "E_M_dom", "I_M_dom",
      # 8-13: number of occupational mosquitoes in each compartment (per individual????)
      "EL_M_occ", "LL_M_occ", "P_M_occ", "S_M_occ", "E_M_occ", "I_M_occ",
      # 14-15: daily EIR due to domestic or occupational mosquitoes
      "EIR_dom", "EIR_occ",
      # 16-21: number of individuals in each compartment
      "S", "I_PCR", "I_LM", "D", "T_", "P",
      # 22: number of people with at least one hypnozoite batch
      "PvHR",
      # 23: number of hypnozoite batches in population
      "PvHR_batch",
      # 24-27: number of incident cases of PCR, LM, symptomatic-untreated or symptomatic-treated disease
      "new_PCR", "new_LM", "new_D", "new_T",
      # 28-30: number of new treatment courses of CQ; number of CQ cases that were relapses or bites
      "new_CQ", "new_CQ_relapse", "new_CQ_bite",
      # 31: number of G6PD tests administered
      "G6PD_tests",
      # 32-34: PQ. number of new treatment courses, number of effective doses administered, number of G6PD-deficient patients treated
      "new_PQ", "PQ_dose_eff", "PQ_G6PD_def",
      # 35-37: TQ. number of new treatment courses, number of effective doses administered, number of G6PD-deficient patients treated
      "new_TQ", "TQ_dose_eff", "TQ_G6PD_def",
      # 38-39: average level of anti-parasite and clinical immunity per individual
      "A_par", "A_clin"
    )
  
  
  # loop over output files -----------------------------------------------------
  
  # initialise lists
  timeline <- list()
  stats    <- list()
  
  for (i in seq_along(list.files(path = current_dir, pattern = "sample_"))) {
    
    out <-
      # read in
      read_tsv(
        file.path(current_dir, paste0("sample_", i, ".tsv")),
        col_names = cols
      ) %>%
      # add sample nr
      mutate(sample_nr = i) %>%
      # transform
      mutate(
        # save time in years
        day  = as.integer(day),
        year = day / 365,
        year_prior = floor(year + 1) %>% as.integer(),
        # other variables to save
        N_pop      = S + I_PCR + I_LM + D + T_ + P,
        prev_PCR   = I_PCR + I_LM + D + T_,
        prev_LM    =         I_LM + D + T_,
        prev_symp  =                D + T_,
        prev_treat =                    T_
      ) %>%
      # drop some columns
      select(-cols[c(2:13, 16:21)])
    
    
    # timeline: smoothen and thin -------------------------------------------- #
    
    timeline[[i]] <-
      out %>%
      # smooth some measurements
      transmute(
        day,
        year,
        year_prior,
        N_pop,
        across(
          .cols  = c(prev_PCR, new_T),
          .fns   = ~ zoo::rollmean(.x, k = 30, align = "center", fill = NA),
          .names = "smooth_{.col}"
        )
      ) %>%
      # thin time points
      filter(day %% 10 == 0) %>%
      # clean
      filter(
        # remove rows where any of the summary measures is NA
        !if_any(starts_with("smooth_"), is.na) & 
          # remove last year of data since the smoothing would not look good
          year_prior < max(year_prior)
      ) %>%
      select(-day, -year_prior)
    
    
    # stats: yearly prevalences (mean) and incidences (sum) ------------------ #
    
    stats[[i]] <-
      out %>%
      # compute cumulative incidence in the years leading up to the timepoints of interest
      group_by(N_pop, sample_nr, year_prior) %>%
      summarise(
        across(
          .cols  = c(prev_PCR),
          .fns   = ~ mean(.x),
          .names = "mean_{.col}"
        ),
        across(
          .cols  = c(new_T, new_CQ_relapse, new_CQ_bite, G6PD_tests, new_CQ, new_PQ, new_TQ),
          .fns   = sum,
          .names = "sum_{.col}"
        )
      ) %>%
      ungroup() %>%
      # clean
      rename(year = year_prior) %>%
      # remove last year of data since it's not complete
      filter(year != max(year))
    
  }
  
  
  # summarise results across runs ----------------------------------------------
  
  # create folder to store partial results
  dir_res <- file.path(dirname(dirname(current_dir)), "results_partial/")
  if (!dir.exists(dir_res))
    dir.create(dir_res)
  
  
  # bind some yearly stats needed to compute effect sizes
  stats %>%
    bind_rows() %>%
    select(N_pop, sample_nr, year, mean_prev_PCR, sum_new_T) %>%
    # add values of varying parameters
    bind_cols(row_varying) %>%
    write_csv(file = paste0(dir_res, "effsize_", basename(current_dir), ".csv"))
  
  
  # summarise remaining yearly stats
  stats %>%
    bind_rows() %>%
    # compute summary stats: median and 95% CrI
    pivot_longer(cols = starts_with(c("sum_", "mean_"))) %>%
    group_by(N_pop, year, name) %>%
    summarise(
      m  = median  (value,                na.rm = TRUE),
      q1 = quantile(value, probs = 0.025, na.rm = TRUE),
      q2 = quantile(value, probs = 0.975, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # compress output
    pivot_wider(
      names_from = name,
      names_sep = "__",
      values_from = c(m, q1, q2)
    ) %>%
    # add values of varying parameters
    bind_cols(row_varying) %>%
    write_csv(file = paste0(dir_res, "stats_", basename(current_dir), ".csv"))
  
  rm(stats)
  
  
  # summarise timelines
  timeline %>%
    bind_rows() %>%
    # compute summary stats: median and 95% CrI
    pivot_longer(cols = starts_with("smooth_")) %>%
    group_by(across(-any_of("value"))) %>%
    summarise(
      m  = median  (value,                na.rm = TRUE),
      q1 = quantile(value, probs = 0.025, na.rm = TRUE),
      q2 = quantile(value, probs = 0.975, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # compress output
    pivot_wider(
      names_from = name,
      names_sep = "__",
      values_from = c(m, q1, q2)
    ) %>%
    # add values of varying parameters
    bind_cols(row_varying) %>%
    write_csv(file = paste0(dir_res, "timeline_", basename(current_dir), ".csv"))
  
  rm(timeline)
  
}