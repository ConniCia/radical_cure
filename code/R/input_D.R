create_input_files_D <- function (step) {
  
  dt_efficacy  <- get_median_efficacy_estimates(taskname = "AON")
  dt_posterior <- get_sample_from_IBM_posterior()
  
  dt <-
    # combine tables of varying and default parameters, and create names for each row
    combine_varying_and_default() %>%
    # combine with drug regimens and remove useless jobs
    combine_drug_regimens() %>%
    # add median hypnozoiticidal efficacy estimates
    bind_cols(dt_efficacy) %>%
    # add sample from IBM posterior
    expand_grid(dt_posterior) %>%
    # transform parameter values
    transform_parameter_values() %>%
    # define input and output file names
    define_file_names(step) %>%
    # create general input files
    create_general_input_files(step) %>%
    # create mosquito input files
    create_mosquito_input_files(step) %>%
    # create intervention input files
    create_intervention_input_files(step)
  
  
  # print parameters sampled from Michael White's IBM posterior
  dt_posterior %>%
    write_tsv(file = file.path(dirX_inp(step), "posterior.tsv"))
  
  
  # print fixed parameters to file
  Filter(function(x) (n_distinct(x) == 1), dt)[1, ] %>%
    write_tsv(file = file.path(dirX_inp(step), "fixed.tsv"))
  
  
  # print varying parameters to file
  Filter(function(x) (n_distinct(x) > 1), dt) %>%
    # remove parameters from posterior
    select(-setdiff(names(dt_posterior), "sample_nr")) %>%
    # move some columns first so bash finds them
    relocate(path_par, path_mosq, path_intv, path_out, sample_nr) %>%
    write_tsv(file = file.path(dirX_inp(step), "varying.tsv"))
  
}


# combine tables of varying and default parameters, and create names for each row
combine_varying_and_default <- function () {
  
  dt <-
    # join default values with varying values: missing values are NA
    bind_rows(dt_default, dt_varying) %>%
    # only keep values that are varying
    lapply(function (col) {
      
      col[which(col == col[1])[-1]] <- NA
      
      return(col)
      
    }) %>%
    bind_cols() %>%
    # remove rows that only contain NAs
    filter(!if_all(everything(), is.na)) %>%
    # remove duplicate rows
    unique() %>%
    # find what pars are varying in each row
    rowwise() %>%
    group_split() %>%
    lapply(function (row) {
      
      # only keep columns that are not NA (these contain non-default parameter values)
      dt <- row %>% select(!where(is.na))
      
      # create strings to use in file naming and creation of sensitivity plots
      if (ncol(dt) == ncol(dt_default)) {
        row <- row %>% mutate(pars_varying = "default", vals_varying = "default")
      } else {
        row <-
          row %>%
          mutate(
            pars_varying = paste(names(dt), collapse = "-"),
            vals_varying = paste(names(dt), dt, sep = "_", collapse = "-")
          )
      }
      
      return(row)
      
    }) %>%
    bind_rows()
  
  
  dt <-
    # replace NA values with default values
    rows_patch(
      dt %>% mutate(rn = row_number()),
      dt_default %>% slice(rep(1, nrow(dt))) %>% mutate(rn = row_number()),
      by = "rn"
    ) %>%
    select(-rn) %>%
    # make strings in pars_varying nicer
    mutate(
      pars_varying =
        pars_varying %>%
        gsub(pattern = "dry_seas-dur_high_seas",     replacement = "seasonality", x = .) %>%
        gsub(pattern = "PQ_7__adhere-PQ_14__adhere", replacement = "PQ__adhere",  x = .)
    )
  
  return(dt)
  
}


# combine with drug regimens and remove useless jobs
combine_drug_regimens <- function (dt) {
  
  dt %>%
    # add 8-AQ drug regimens
    expand_grid(drug_regimen = label_drug_regimen %>% filter(! is.na(colour)) %>% pull(name)) %>%
    # remove "useless" jobs:
    # PQ adherence is not used in No_8AQs, AQ_ideal and TQ regimens
    filter(! ( grepl("AQ|TQ", drug_regimen) & grepl("PQ__adhere", pars_varying) ) ) %>%
    # G6PD prevalence and CYP2D6 prevalence are not used in No_8AQs and AQ_ideal regimens
    filter(! ( grepl("AQ",    drug_regimen) & grepl("G6PD_prev|CYP2D6_prev", pars_varying) ) )
  
}


# get median hypnozoiticidal efficacy estimates
get_median_efficacy_estimates <- function (taskname) {
  
  c("A", "B", "C") %>%
    lapply(function (step) {
      
      read_estimates(step, "MCMC", taskname) %>%
        handle_burnin(action = "remove") %>%
        # select efficacy columns
        select(starts_with("eff")) %>%
        # compute median
        summarise(across(everything(), median))
      
    }) %>%
    bind_cols()
  
}


# get sample of the posterior parameter draws obtained from the calibration of 
# the PvIBM by [White et al](http://dx.doi.org/10.1038/s41467-018-05860-8)
get_sample_from_IBM_posterior <- function () {
  
  set.seed(1)
  
  read_csv(file.path(dir_raw, "PvMod_posterior_draws.csv")) %>%
    # sample 'n_sample' parameter combinations from posterior
    slice_sample(n = n_sample) %>%
    # define sample_nr
    mutate(sample_nr = row_number())
  
}


# transform parameter values
transform_parameter_values <- function (dt) {
  
  dt %>%
    mutate(
      # rename (keeping old and new columns)
      CM_CQ_eff     = CQ__eff,
      CM_CQ_eff_wPQ = CQ__eff_wPQ,
      CM_CQ_proph   = CQ__proph,
      A_PCR_50pc    = A_d_PCR_50pc,
      K_PCR         = K_d_PCR,
      P_mat         = P_MI,
      d_mat         = d_MI,
      # transform durations into rates
      ff      = 1 / d_ff,
      gamma_L = 1 / d_gamma,
      r_LM    = 1 / d_LM,
      # transform yearly into daily quantities
      EIR_equil = yEIR / 365,
    ) %>%
    # compute kappa_seas
    group_by(dur_high_seas) %>%
    mutate(kappa_seas = find_kappa(dur_high_seas[1])) %>%
    ungroup() %>%
    # limit number of decimal digits
    mutate(across(where(is.numeric), ~ round(.x, digits = 9)))
  
}


# define input and output file names
define_file_names <- function (dt, step) {
  
  dt %>%
    group_split(vals_varying) %>%
    lapply(function (x) {
      
      vv <- x$vals_varying[1]
      pv <- x$pars_varying[1]
      
      if (pv == "CM_cover") {
        name_par  <- vv
        name_mosq <- "default"
        name_intv <- vv
      } else if (pv == "PQ__adhere") {
        name_par  <- "default"
        name_mosq <- "default"
        name_intv <- vv
      } else if (pv == "seasonality") {
        name_par  <- "default"
        name_mosq <- vv
        name_intv <- "default"
      } else if (pv == "yEIR-PQ__adhere") {
        name_par  <- gsub(x = vv,  pattern = "(-).*", replacement = "")
        name_mosq <- "default"
        name_intv <- gsub(x = vv,  pattern = "(yEIR)[^-]*-", replacement = "")
      } else if (pv == "yEIR-seasonality") {
        name_par  <- gsub(x = vv,  pattern = "(-).*", replacement = "")
        name_mosq <- gsub(x = vv,  pattern = "(yEIR)[^-]*-", replacement = "")
        name_intv <- "default"
      } else {
        name_par  <- vv
        name_mosq <- "default"
        name_intv <- "default"
      }
      
      x %>%
        # create file paths
        mutate(
          vals_varying = NULL,
          path_par  = file.path(dirX_inp(step), "par",  paste0(name_par, "-sample_", sample_nr, ".txt")),
          path_mosq = file.path(dirX_inp(step), "mosq", paste0(name_mosq, ".txt")),
          path_intv = file.path(dirX_inp(step), "intv", paste0(drug_regimen, "-", name_intv, ".txt")),
          path_out  = file.path(dirX_out(step), paste0(drug_regimen, "-", vv), paste0("sample_", sample_nr, ".tsv"))
        )
      
    }) %>%
    bind_rows()
  
}


# create general input files
create_general_input_files <- function (dt, step) {
  
  file_type <- "par"
  
  # create input folder
  dir_inp <- dirname(dt$path_par[1])
  if (!dir.exists(dir_inp))
    dir.create(dir_inp, recursive = TRUE)
  
  
  dt %>%
    # nest
    nest(.by = c(get_empty_table(file_type)$name, path_par)) %>%
    # create files
    rowwise() %>%
    group_walk(.f = function (row, ...) {
      
      # fill empty input table
      left_join(
        row %>% select(where(is.numeric)) %>% pivot_longer(cols = everything()),
        get_empty_table(file_type),
        by = "name"
      ) %>%
        # print to file
        write_tsv(file = row$path_par, col_names = FALSE)
      
    }) %>%
    # unnest
    unnest(cols = data)
  
}


# create mosquito input files
create_mosquito_input_files <- function (dt, step) {
  
  file_type <- "mosq"
  
  # create input folder
  dir_inp <- dirname(dt$path_mosq[1])
  if (!dir.exists(dir_inp))
    dir.create(dir_inp, recursive = TRUE)
  
  
  dt %>%
    # nest
    nest(.by = c(get_empty_table(file_type)$name, path_mosq)) %>%
    # create files
    rowwise() %>%
    group_walk(.f = function (row, ...) {
      
      # fill empty input table
      left_join(
        row %>% select(where(is.numeric)) %>% pivot_longer(cols = everything()),
        get_empty_table(file_type),
        by = "name"
      ) %>%
        # print to file
        write_tsv(file = row$path_mosq, col_names = FALSE)
      
    }) %>%
    # unnest
    unnest(cols = data)
  
}


# create intervention input files
create_intervention_input_files <- function (dt, step) {
  
  file_type <- "intv"
  
  # create input folder
  dir_inp <- dirname(dt$path_intv[1])
  if (!dir.exists(dir_inp))
    dir.create(dir_inp, recursive = TRUE)
  
  
  dt %>%
    # create input file using a subset of grouped columns
    nest(.by = c(drug_regimen, CM_cover, intervention_time, contains("__"), path_intv)) %>%
    rowwise() %>%
    group_walk(.f = function (row, ...) {
      
      # prep intervention data
      if (row$drug_regimen == "No_8AQs") {
        tbl <-
          tibble(
            years     = row$intervention_time,
            # CQ
            CM0_cover = row$CM_cover,
            CM0_eff   = row$CQ__eff,
            CM0_proph = row$CQ__proph
          )
      } else if (row$drug_regimen == "AQ_ideal") {
        tbl <-
          tibble(
            years              = row$intervention_time,
            CM1_cover          = row$CM_cover,
            # CQ
            CM1_CQ_eff         = row$CQ__eff,
            CM1_CQ_eff_wPQ     = row$CQ__eff_wPQ,
            CM1_CQ_proph       = row$CQ__proph,
            # AQ_ideal
            CM1_PQ_eff         = row$AQ_ideal__eff,
            CM1_PQ_proph       = row$AQ_ideal__proph,
            CM1_PQ_adhere      = row$AQ_ideal__adhere,
            CM1_PQ_lowage      = row$AQ_ideal__lowage,
            CM1_PQ_G6PD_risk   = row$AQ_ideal__G6PD_risk,
            CM1_PQ_CYP2D6_risk = row$AQ_ideal__CYP2D6_risk,
            CM1_PQ_preg_risk   = row$AQ_ideal__preg_risk,
            # additional G6PD screening
            CM1_G6PD_test      = row$AQ_ideal__G6PD_test
          )
      } else if (grepl("PQ_", row$drug_regimen)) {
        tbl <-
          tibble(
            years              = row$intervention_time,
            CM1_cover          = row$CM_cover,
            # CQ
            CM1_CQ_eff         = row$CQ__eff,
            CM1_CQ_eff_wPQ     = row$CQ__eff_wPQ,
            CM1_CQ_proph       = row$CQ__proph,
            # PQ
            CM1_PQ_eff         = 
              case_match(
                row$drug_regimen,
                "PQ_highdose7"  ~ row$eff__PQ_highdose7,
                "PQ_highdose14" ~ row$eff__PQ_highdose14,
                .default        = row$eff__PQ_lowdose
              ),
            CM1_PQ_proph       = 
              ifelse(
                endsWith(row$drug_regimen, "7"),
                row$PQ_7__proph,
                row$PQ_14__proph
              ),
            CM1_PQ_adhere      =
              ifelse(
                endsWith(row$drug_regimen, "7"),
                row$PQ_7__adhere,
                row$PQ_14__adhere
              ),
            CM1_PQ_lowage      = row$PQ__lowage,
            CM1_PQ_G6PD_risk   = row$PQ__G6PD_risk,
            CM1_PQ_CYP2D6_risk = row$PQ__CYP2D6_risk,
            CM1_PQ_preg_risk   = row$PQ__preg_risk,
            # additional G6PD screening
            CM1_G6PD_test      = row$PQ__G6PD_test
          )
      } else if (grepl("PQrestr_", row$drug_regimen)) {
        # restricted high-dose PQ over 7 days: assume it has the same G6PD restrictions as TQ
        tbl <-
          tibble(
            years              = row$intervention_time,
            CM2_cover          = row$CM_cover,
            # CQ
            CM2_CQ_eff         = row$CQ__eff,
            CM2_CQ_eff_wPQ     = row$CQ__eff_wPQ,
            CM2_CQ_proph       = row$CQ__proph,
            # PQ_lowdose7
            CM2_PQ_eff         = row$eff__PQ_lowdose,
            CM2_PQ_proph       = row$PQ_7__proph,
            CM2_PQ_adhere      = row$PQ_7__adhere,
            CM2_PQ_lowage      = row$PQ__lowage,
            CM2_PQ_G6PD_risk   = row$PQ__G6PD_risk,
            CM2_PQ_CYP2D6_risk = row$PQ__CYP2D6_risk,
            CM2_PQ_preg_risk   = row$PQ__preg_risk,
            # PQrestr_highdose7
            CM2_TQ_eff         = row$eff__PQ_highdose7,
            CM2_TQ_proph       = row$PQ_7__proph,
            CM2_TQ_adhere      = row$PQ_7__adhere,
            CM2_TQ_lowage      = row$PQ__lowage,
            CM2_TQ_G6PD_risk   = row$PQ__G6PD_risk,
            CM2_TQ_CYP2D6_risk = row$PQ__CYP2D6_risk,
            CM2_TQ_preg_risk   = row$PQ__preg_risk,
            # additional G6PD screening
            CM2_G6PD_test      = row$PQ__G6PD_test
          )
      } else if (grepl("TQ_", row$drug_regimen)) {
        tbl <-
          tibble(
            years              = row$intervention_time,
            CM2_cover          = row$CM_cover,
            # CQ
            CM2_CQ_eff         = row$CQ__eff,
            CM2_CQ_eff_wPQ     = row$CQ__eff_wPQ,
            CM2_CQ_proph       = row$CQ__proph,
            # PQ_lowdose7
            CM2_PQ_eff         = row$eff__PQ_lowdose,
            CM2_PQ_proph       = row$PQ_7__proph,
            CM2_PQ_adhere      = row$PQ_7__adhere,
            CM2_PQ_lowage      = row$PQ__lowage,
            CM2_PQ_G6PD_risk   = row$PQ__G6PD_risk,
            CM2_PQ_CYP2D6_risk = row$PQ__CYP2D6_risk,
            CM2_PQ_preg_risk   = row$PQ__preg_risk,
            # TQ
            CM2_TQ_eff         = 
              ifelse(
                row$drug_regimen == "TQ_lowdose",
                row$eff__TQ_lowdose,
                row$eff__TQ_highdose
              ),
            CM2_TQ_proph       = 
              ifelse(
                row$drug_regimen == "TQ_lowdose",
                row$TQ_lowdose__proph,
                row$TQ_highdose__proph
              ),
            CM2_TQ_adhere      = row$TQ__adhere,
            CM2_TQ_lowage      = row$TQ__lowage,
            CM2_TQ_G6PD_risk   = row$TQ__G6PD_risk,
            CM2_TQ_CYP2D6_risk = row$TQ__CYP2D6_risk,
            CM2_TQ_preg_risk   = row$TQ__preg_risk,
            # additional G6PD screening
            CM2_G6PD_test      = row$TQ__G6PD_test
          )
      }
      
      
      # fill empty input table
      get_empty_table(file_type) %>%
        mutate(temp_col = -10) %>%
        # add input values
        pivot_wider(values_from = temp_col) %>%
        bind_rows(purrr::reduce(list(tbl), dplyr::full_join, by = "years")) %>%
        # arrange by years
        arrange(years) %>%
        # replace NA with -1
        replace(is.na(.), -1) %>%
        # transform back to tibble
        t() %>%
        data.frame() %>%
        tibble::rownames_to_column() %>%
        # clean up
        select(-X1) %>%
        mutate(last_col = -1) %>%
        # print to file
        write_delim(file = row$path_intv, col_names = FALSE)
      
    }) %>%
    unnest(cols = data)
  
}

