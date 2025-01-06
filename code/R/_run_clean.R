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
current_dir     <- args[2]


# source functions
setwd(path_on_cluster)
source("code/R/_header.R")


clean_D(current_dir)
