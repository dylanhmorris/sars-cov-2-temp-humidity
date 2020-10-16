#!/usr/bin/env Rscript

###################################
## filename: literature_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of titer decay
## that uses data taken from the literature
####################################

hyperparam_list <- list(
    log_hl_prior_mean = -2,
    log_hl_prior_sd = 4,
    mode_sd_sample = 0.6,
    sd_sd_sample = 0.2,
    debug = FALSE)

fixed_seed <- 23032
inits_seed <- 32
