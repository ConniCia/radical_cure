#
# ATTENTION
#
# This R script is called by bash and run directly on the cluster!
#
# Don't source or run this locally, it will take way too much time!
#


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

path_on_cluster <- args[1]
step            <- args[2]
input_file      <- args[3]
N_LHS    <- as.integer(args[4])
N_MCMC   <- as.integer(args[5])
n_chains <- as.integer(args[6])
burnin   <- as.numeric(args[7])


# source functions
setwd(path_on_cluster)
source("code/R/_header.R")


# run fit
if (step == "A") {
  fit_A(step, input_file, N_LHS, N_MCMC)
} else {
  fit_BC(step, input_file, N_LHS, N_MCMC)
}
