# Load packages
library(dplyr)      # all of the code!
library(tidyr)      # all of the code!
library(readr)      # all of the code!
library(purrr)      # all of the code! [map()]
library(parallel)   # [mclapply(), detectCores()] run lapply in parallel
library(lhs)        # [randomLHS()] sample of points for Latin Hypercube Sampling
library(ggplot2)    # plots
library(zoo)        # clean results


require(odin)       # ODE system
require(stringr)    # [str_remove()]
require(mcmcse)     # [ess(), multiESS(), minESS()] compute effective sample size
require(grid)       # plots
require(gridExtra)  # plots
require(scales)     # numeric axis labels
require(latex2exp)  # [TeX()] plot mathematical symbols in figures
require(kableExtra) # nice table styling
require(forcats)    # [fct_rev()]
require(ggpattern)  # [geom_col_pattern()]
# require(plyr)       # [revalue()] -> problematic to load plyr AFTER dyplr!!!!


# Console options
mc.cores <- parallel::detectCores() - 1L # for use in mclapply()
options(
  scipen = 999,                    # Disable scientific notation in format()
  digits = 7,                      # Nr of significant digits used in format()
  dplyr.summarise.inform = FALSE,  # Disable dplyr messages
  readr.show_col_types   = FALSE,  # Disable message when reading in csv
  readr.show_progress    = FALSE,  # Avoid empty lines below call to read_csv()
  odin.verbose           = FALSE   # Disable message when compiling odin models
)


# Source all R scripts that don't start with an underscore
list.files(path = "code/R", full.names = TRUE) %>%
  grep(pattern = "R/_", value = TRUE, invert = TRUE) %>%
  sapply(source) %>%
  invisible()


# Define paths
address_cluster <- "maestro.pasteur.fr"
path_on_cluster <- "/pasteur/zeus/projets/p02/Projets_Mph_HPC/Conni/RP13/"


# Define directories
dir_raw  <- "data/raw"
dir_diag <- "data/diagnostics"
dir_fig  <- "data/figures"
dir_tab  <- "data/tables"

steps <- LETTERS[1:4]

dirX     <- function (step) file.path("data", step)
dirX_inp <- function (step) file.path("data", step, "input")
dirX_out <- function (step) file.path("data", step, "output")
dirX_res <- function (step) file.path("data", step, "results")

create_dirs <- function (steps) {
  
  suppressWarnings(sapply(
    c(dirX_inp(steps), dirX_out(last(steps)), dirX_res(steps), dir_diag, dir_fig, dir_tab),
    dir.create,
    recursive = TRUE
  ))
  
}

create_dirs(steps)


# Global figure theme
theme_global <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                          base_rect_size = base_size/22) {
  
  half_line <- base_size/2
  
  theme_light(base_size, base_family, base_line_size, base_rect_size) %+replace%
    theme(
      panel.border = element_rect(fill = NA, colour = NA),
      strip.background = element_rect(fill = "lightgrey", colour = NA),
      strip.text = element_text(colour = "black", size = rel(0.8), margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)), 
      axis.ticks = element_line(colour = "black"),
      legend.margin = margin(0,0,0,0),
      complete = TRUE
    )
  
}

