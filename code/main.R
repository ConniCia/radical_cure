############################################################################## #
##                                                                            ##
##        EVALUATE THE EFFECT OF INTRODUCING PRIMAQUINE OR TAFENOQUINE        ##
##                        IN P. VIVAX CASE MANAGEMENT                         ##
##                                                                            ##
##                Authors: Constanze Ciavarella, Michael White                ##
##      This code: Constanze Ciavarella (constanze.ciavarella@pasteur.fr)     ##
##                                                                            ##
##                  (for details, please consult the README)                  ##
##                                                                            ##
############################################################################## #


# ---------------------------------------------------------------------------- #
#  Preliminaries                                                               #
# ---------------------------------------------------------------------------- #

## Set up workspace
source("code/R/_header.R")

## Global options
burnin   <- 0.25  # what proportion of MCMC chains to throw away as burnin
n_sample <- 100L  # nr of parameter combinations to sample from posteriors
n_chains <- 6     # nr of MCMC chains

## Default parameters (tibble)
dt_default <- default_parameters()

## Varying parameters (tibble)
dt_varying <- varying_parameters()


# ---------------------------------------------------------------------------- #
#  Run cluster jobs                                                            #
# ---------------------------------------------------------------------------- #

## A: Run MCMC to fit efficacy of PQ_highdose7 and PQ_highdose14
run(step = "A", N_LHS = 100000, N_MCMC = 50000)

## B: Run MCMC to fit efficacy of PQ_lowdose (7 and 14 days confounded)
run(step = "B", N_LHS = 10000, N_MCMC = 50000)

## C: Run MCMC to fit efficacy of TQ_lowdose and TQ_highdose
run(step = "C", N_LHS = 10000, N_MCMC = 50000)

## D: Run IBM to simulate introduction of 8AQ in case management
run(step = "D")

## D: Clean IBM output
run_cleaning(step = "D")

## Copy latest cluster results
copy_cluster_results()


monitor_jobs(compact = TRUE)
monitor_jobs(compact = FALSE)
cancel_jobs()


# ---------------------------------------------------------------------------- #
#  Create figures and tables                                                   #
# ---------------------------------------------------------------------------- #

## graphic devices
device_SI   <- "png"
device_main <- c(device_SI, "pdf")

## A, B, C: diagnostic plots
create_diagnostic_plots()

## A, B, C: paper figures
create_figures_ABC()

## D: paper figures
create_figures_D()

## all steps: paper tables
create_tables()
