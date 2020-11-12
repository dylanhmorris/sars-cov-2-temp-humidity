#!/usr/bin/env Rscript

########################################
## filename: macros-energy-comparison.R
## author: Dylan Morris <dhmorris@princeton.edu>
## description: create macros for the posterior
## of the difference between E_a_sol and E_a_eff
#######################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for spread_draws()
    'virusenv'    # project functions
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

#################################
## read command line args
#################################
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
dat_path <- args[1]
evap_dat_path <- args[2]
meas_results_path <- args[3]
mod_results_path <- args[4]
outpath <- args[5]

#################################
# read in needed data
#################################

dat <- read_data_for_plotting(dat_path)

meas_results_chains <- readRDS(meas_results_path)
mod_results_chains <- readRDS(mod_results_path)


parse_chains <- function(chains){
    result <- chains %>%
        spread_draws(log_k_30_20_ratio_dry,
                     log_k_30_20_ratio_wet,
                     E_a_dry,
                     E_a_wet) %>%
        mutate(log_ratio_ratio =
                   log_k_30_20_ratio_dry - log_k_30_20_ratio_wet,
               E_a_diff = E_a_dry - E_a_wet,
               E_a_percent = E_a_diff / E_a_dry)
    
    return(result)
}

meas_results_draws <- meas_results_chains %>%
    parse_chains() %>%
    mutate(model = "measured concentration")

mod_results_draws <- mod_results_chains %>%
    parse_chains() %>%
    mutate(model = "modeled concentration")


#################################
## macro setup
#################################
macro_list <- list()
macro_template <- "\\newcommand{\\%s}{%s}\n"

    cat('saving macros to', outpath, '...\n')

macro_list["MedianEaDiffMeasured"] <-
    100 * median(meas_results_draws$E_a_percent)
macro_list["QLowEaDiffMeasured"] <-
    100 * quantile(meas_results_draws$E_a_percent, 0.025)
macro_list["QHighEaDiffMeasured"] <-
    100 * quantile(meas_results_draws$E_a_percent, 0.975)

macro_list["MedianEaDiffModeled"] <-
    100 * median(mod_results_draws$E_a_percent)
macro_list["QLowEaDiffModeled"] <-
    100 * quantile(mod_results_draws$E_a_percent, 0.025)
macro_list["QHighEaDiffModeled"] <-
    100 * quantile(mod_results_draws$E_a_percent, 0.975)

if (file.exists(outpath))
    file.remove(outpath)

for(param_name in names(macro_list)){
    cat(sprintf(macro_template,
                param_name,
                macro_list[param_name]),
        file = outpath,
        append = TRUE)
}


warnings()
