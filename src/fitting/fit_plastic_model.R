#!/usr/bin/env Rscript

###################################
## filename: fit_plastic_model.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: fit the specified stan model
## for virus stability on plastic
## as a function of temperature and
## humidity
####################################

script_packages <- c(
    'rstan',     # stan interface
    'parallel',  # parallelized MCMC
    'readr',     # csv read-in
    'magrittr',  # for pipe operator %>%
    'dplyr',     # for filter()
    'tidyr'      # for drop_na()
)


## set up hyperparameters for models

debug <- FALSE ## set to TRUE to diagnose sampler problems
strict <- FALSE ## set to TRUE to reject divergent transition

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <- args[1]
data_path <- args[2]
hyperparam_path <- args[3]
mcmc_output_path <- args[4]
prior_check <- (as.logical(args[5]) &
                (!is.na(as.logical(args[5]))))
empirical_salt <- (as.logical(args[6]) &
                   (!is.na(as.logical(args[6]))))




## read data
cat('reading in data from file ', data_path, ' ...\n')
dat <- read_csv(data_path,
                col_types = cols())
cat('data loaded successfully!\n')

cat('reading in hyperparameters from file ', hyperparam_path, ' ...\n')
source(hyperparam_path)
cat('hyperparameters loaded successfully!\n')

## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 2000
nchains <- n_cores
adapt_d <- 0.95
max_tree <- 10
set.seed(inits_seed) ## R's rng set for random inits

if (debug) {
    nchains <- 1
    niter <- 1000
}

## calculate munged data
exp_dat <- dat %>%
    distinct(experiment_id, .keep_all = TRUE)

if(!"equilibrium_mass" %in% names(exp_dat)){
    exp_dat$equilibrium_mass <- NA
    exp_dat$initial_mass <- NA
}

titer_dat <- dat %>%
    distinct(titer_id,
             .keep_all = TRUE)


## defined values
n_titers <-  max(titer_dat$titer_id)
n_experiments <- max(exp_dat$experiment_id)
ref_exp_id <- exp_dat %>%
    filter(temperature == 22 &
           humidity == 65) %>%
    select(experiment_id) %>%
    first()
n_total_datapoints <- length(dat$virus_detect)

## null values for prior checks
n_used_datapoints <- 0
n_used_titers <- 0

## non null values for fitting
if(!prior_check){
    n_used_titers <- n_titers
    n_used_datapoints <- n_total_datapoints
}

observation_data_list <- list(
    n_total_datapoints = n_total_datapoints,
    n_used_datapoints = n_used_datapoints,
    n_experiments =  n_experiments,
    n_titers = n_titers,
    n_used_titers = n_used_titers,
    n_transient_classes = max(exp_dat$transient_class_id),
    titer_id = dat$titer_id,
    dilution = dat$dilution,
    well_status = dat$virus_detect,
    experiment_id = dat$experiment_id,
    temperature = exp_dat$temperature,
    transient_class_id = exp_dat$transient_class_id,
    rel_humidity = exp_dat$humidity / 100,
    titer_experiment_id = titer_dat$experiment_id,
    titer_times = titer_dat$time,
    equilibrium_mass = exp_dat$equilibrium_mass,
    initial_mass = exp_dat$initial_mass,
    ref_exp_id = ref_exp_id,
    empirical_salt = empirical_salt)

init_val <- "random"

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

cat('\nsaving mcmc samples to', mcmc_output_path, '\n')

saveRDS(fit, mcmc_output_path)

warnings()
