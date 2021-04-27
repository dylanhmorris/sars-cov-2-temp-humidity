#!/usr/bin/env Rscript

#####################################
## name: chain_diagnostics.R
## author: Dylan Morris
##
## read in mcmc chains
## and output diagnostics
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot)) # for rhat, neff_ratio

## read command line args
args <- commandArgs(trailingOnly = TRUE)
n_models <- length(args) - 1
mcmc_fit_paths <- args[1:n_models]
outpath <- args[n_models + 1]

rhat_max <- function(chains){
    rhat_vals <- rhat(chains)

    if (any(is.nan(rhat_vals))) {
        warning(paste0(
            "Warning: removing NaN rhat values. This may be benign ",
            "if, for example, there is a parameter that is not sampled ",
            "because it is fully constrained"))
        rhat_vals <- rhat_vals[!is.nan(rhat_vals)]
    }

    rhat_max_parm <- names(which.max(rhat_vals))
    rhat_max_val <- max(rhat_vals)

    return (list(rhat_max_parm, rhat_max_val))

}

neff_min <- function(chains){

    neff_ratios <- neff_ratio(chains)
    
    if (any(is.nan(neff_ratios))) {
        warning(paste0(
            "Warning: removing NaN neff_ratio ",
            "values. This may be benign ",
            "if, for example, there is a parameter ",
            "that is not sampled ",
            "because it is fully constrained"))
        neff_ratios <- neff_ratios[!is.nan(neff_ratios)]
    }

    neff_min_parm <- names(which.min(neff_ratios))
    neff_min_val <- min(neff_ratios)
    return (list(neff_min_parm, neff_min_val))

}

## create empty structure
models <- rep('', n_models)
rhat_max_parms <- rep('', n_models)
rhat_max_vals <- rep(NA, n_models)
neff_min_parms <- rep('', n_models)
neff_min_vals <- rep(NA, n_models)
n_divs <- rep(NA, n_models)

## get mcmc chains
for (k in 1:length(mcmc_fit_paths)) {

    path <- mcmc_fit_paths[k]
    model_name <- path
    
    cat(sprintf("\nDiagnostics for model %s...\n",
                model_name))
    fit <- readRDS(path)
    chains <- fit
    rhat_maxes <- rhat_max(chains)
    neff_mins <- neff_min(chains)

    models[k] <- model_name
    rhat_max_parms[k] <- rhat_maxes[[1]]
    rhat_max_vals[k] <- rhat_maxes[[2]]
    neff_min_parms[k] <- neff_mins[[1]]
    neff_min_vals[k] <- neff_mins[[2]]

    sampler_params <- get_sampler_params(fit,
                                         inc_warmup = FALSE)
    divs <- sum(sapply(sampler_params,
                       function(x) sum(x[, "divergent__"])))

    n_divs[k] <- divs
}

diagnostic_data <- data.frame(
    models = models,
    rhat_max_parm = rhat_max_parms,
    rhat_max_val = rhat_max_vals,
    neff_min_parm = neff_min_parms,
    neff_min_val = neff_min_vals,
    n_divergent = n_divs)

write.csv(diagnostic_data, outpath)
