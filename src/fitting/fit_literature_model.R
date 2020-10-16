#!/usr/bin/env Rscript

#####################################
## name: fit_literature_model.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## fit models of inactivation of
## for virus data from the literature
####################################


suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(dplyr))

## functions

fit_model <- function(model_src_path,
                      data_path,
                      hyperparam_path,
                      output_path,
                      prior_check,
                      strict = FALSE,
                      debug = FALSE) {

    cat('reading in data from file ',
        data_path, ' ...\n')
    dat <- read_csv(data_path,
                    col_types = cols())
    cat('data loaded successfully!\n')

    cat('reading in hyperparameters from file ',
        hyperparam_path, ' ...\n')
    source(hyperparam_path)
    cat('hyperparameters loaded successfully!\n')

    ## set stan options
    n_cores <- parallel::detectCores()
    options(mc.cores = n_cores)
    rstan_options(auto_write = TRUE)
    niter <- 3000
    nchains <- n_cores
    adapt_d <- 0.95
    max_tree <- 10
    set.seed(inits_seed) ## R's rng set for random inits
    init_val <- "random"

    if (debug) {
        nchains <- 1
        niter <- 1000
    }

    ## calculate munged data
    exp_dat <- dat %>%
        distinct(experiment_id, .keep_all = TRUE)

    ## null values for prior checks
    n_measurements_used <- 0

    ## non null values for fitting
    if(!prior_check){
        n_measurements_used <- max(dat$measurement_id)
    }

    observation_data_list <- list(
        n_measurements = max(dat$measurement_id),
        n_measurements_used = n_measurements_used,
        n_experiments =  max(exp_dat$experiment_id),
        n_experiment_classes =  max(dat$experiment_class_id),
        experiment_id = dat$experiment_id,
        experiment_class_id = dat$experiment_class_id,
        log10_fraction_measured = dat$log10_fraction_viable,
        time = dat$elapsed_time,
        log10_fraction_LOD = dat$log10_LOD_fraction)


    ###############################
    ## Compile, fit, and save model
    ###############################
    ## pass stan data and hyperparams
    stan_data <- c(
        observation_data_list,
        hyperparam_list)

    fit <- stan(
        model_src_path,
        data = stan_data,
        iter = niter,
        seed = fixed_seed,
        chains = nchains,
        init = init_val,
        control = list(max_treedepth = max_tree,
                       adapt_delta = adapt_d))


    ## check that sampled correctly and only save
    ## if so
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)

    if(is.null(sampler_params)){
        stop("Stan model failed to sample")
    }

    divs <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

    if(strict){
        if(divs > 0){
            stop("Divergent transitions when sampling. Check priors and parametrization")
        }
    }

    cat('\nsaving mcmc samples to', output_path, '\n')

    saveRDS(fit, output_path)

}


## run script
args <- commandArgs(trailingOnly = TRUE)
model_src_path <- args[1]
data_path <- args[2]
hyperparam_path <- args[3]
mcmc_output_path <- args[4]
prior_check <- (as.logical(args[5]) &
                (!is.na(as.logical(args[5]))))

fit_model(model_src_path,
          data_path,
          hyperparam_path,
          mcmc_output_path,
          prior_check)

warnings()
