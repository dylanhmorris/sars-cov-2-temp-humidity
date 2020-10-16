#!/usr/bin/env Rscript

###################################
## filename: infer_titer_hypers.R
 ## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of titer inference
####################################

virology_hypers <- list(
    titer_prior_mean = 2.5,
    titer_prior_sd = 4)

flags <- list(
    debug = FALSE)

hyperparam_list <- c(
    virology_hypers,
    flags)

fixed_seed <- 76543
inits_seed <- 1632
