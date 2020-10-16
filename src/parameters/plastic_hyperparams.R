#!/usr/bin/env Rscript

###################################
## filename: plastic_hyperparams.R
 ## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of titer decay
## that uses raw well data and fits
## a dependence on temperature and humidity
####################################

virology_hypers <- list(
    intercept_prior_mean = 2.5, # official dose was 10^5 TCID50
    intercept_prior_sd = 1,
    mode_sd_intercept = 0,
    sd_sd_intercept = 0.5)

inactivation_hypers <- list(
    log_hl_20_prior_mean = log(24),
    log_hl_20_prior_sd = 1.25,
    log_hl_prior_mean = log(6),
    log_hl_prior_sd = 2,
    log_hl_transient_prior_mean = log(24),
    log_hl_transient_prior_sd = 1.25,
    log_k_30_20_ratio_prior_mean = 0,
    log_k_30_20_ratio_prior_sd = 1)

evaporation_hypers <- list(
    ERH = 0.45,
    mean_drying_time = 10,
    sd_drying_time = 10,
    mode_sd_evap = 0,
    sd_sd_evap = 1,
    mode_sd_logconc = 0,
    sd_sd_logconc = 1,
    log_initial_mass_solute_per_thousand_prior_mean = log(11),
    log_initial_mass_solute_per_thousand_prior_sd = 0.33,
    mean_log_shape_parameter_concavity = 0,
    sd_log_shape_parameter_concavity = 0.33,
    mean_log_shape_parameter_steepness = 0,
    sd_log_shape_parameter_steepness = 0.33)


flags <- list(
    functional_form = 1,
    debug = FALSE)

hyperparam_list <- c(
    virology_hypers,
    inactivation_hypers,
    evaporation_hypers,
    flags)

fixed_seed <- 923
inits_seed <- 1112
