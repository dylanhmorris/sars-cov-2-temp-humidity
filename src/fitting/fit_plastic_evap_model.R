#!/usr/bin/env Rscript

###################################
## filename: fit_plastic_evap_model.R
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
args <- commandArgs(trailingOnly = TRUE)
model_src_path <- args[1]
data_path <- args[2]
evap_data_path <- args[3]
hyperparam_path <- args[4]
mcmc_output_path <- args[5]
prior_check <- (as.logical(args[6]) &
                (!is.na(as.logical(args[6]))))
empirical_salt <- (as.logical(args[7]) &
                   (!is.na(as.logical(args[7]))))
shared_E_a <- (as.logical(args[8]) &
               (!is.na(as.logical(args[8]))))




## read data
cat('reading in data from file ', data_path, ' ...\n')
dat <- read_csv(data_path,
                col_types = cols())

cat('reading in data from file ', evap_data_path, ' ...\n')
evap_dat <- read_csv(evap_data_path,
                     col_types = cols())
print(evap_dat)
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
adapt_d <- 0.9

if(grepl("sars-mers", mcmc_output_path)){
    adapt_d <- 0.99
}
   
max_tree <- 10
set.seed(inits_seed) ## R's rng set for random inits

if (debug) {
    nchains <- 1
    niter <- 1000
}

## calculate munged data

titer_dat <- dat %>%
    distinct(titer_id,
             .keep_all = TRUE)


exp_dat <- dat %>%
    distinct(experiment_id, .keep_all = TRUE)

evap_dat <- evap_dat %>%
    filter(measured_mass > equilibrium_mass) %>%
    inner_join(exp_dat %>% select(temperature,
                                  humidity,
                                  evap_class_id),
               by = c("temperature", "humidity")) %>%
    ungroup()

evap_class_dat <- evap_dat %>%
    distinct(evap_class_id, .keep_all = TRUE)

exp_dat <- exp_dat %>%
    inner_join(evap_class_dat %>%
               select(temperature,
                      humidity,
                      experiment_evap_class = evap_class_id),
               by = c("temperature", "humidity"))

print(evap_dat %>% distinct(temperature, humidity, evap_class_id))
print(evap_class_dat %>% select(temperature, humidity, evap_class_id, initial_mass, equilibrium_mass))
print(exp_dat %>% select(temperature, humidity, experiment_id, experiment_evap_class))


## defined values
n_total_datapoints <- length(dat$virus_detect)
n_titers <-  max(titer_dat$titer_id)
n_experiments <- max(exp_dat$experiment_id)
n_evap_datapoints <- length(evap_dat$measured_mass)

## default use all data
n_used_titers <- n_titers
n_used_datapoints <- n_total_datapoints
n_used_evap_datapoints <- n_evap_datapoints
n_evap_classes <- max(evap_dat$evap_class_id)
n_transient_classes <- max(exp_dat$transient_class_id)

## null values for prior checks
if(prior_check){
    n_used_titers <- 0
    n_used_datapoints <- 0
    n_used_evap_datapoints <- 0
}

## compose data lists
observation_data_list <- list(
    n_total_datapoints = n_total_datapoints,
    n_used_datapoints = n_used_datapoints,
    n_experiments =  n_experiments,
    n_transient_classes = n_transient_classes,
    n_titers = n_titers,
    n_used_titers = n_used_titers,
    titer_id = dat$titer_id,
    dilution = dat$dilution,
    well_status = dat$virus_detect,
    experiment_id = dat$experiment_id)

experiment_data_list <- list(
    temperature = exp_dat$temperature,
    rel_humidity = exp_dat$humidity / 100,
    transient_class_id = exp_dat$transient_class_id,
    experiment_evap_class = exp_dat$experiment_evap_class)

titer_data_list <- list(
    titer_experiment_id = titer_dat$experiment_id,
    titer_times = titer_dat$time)

evap_data_list <- list(
    n_evap_datapoints = n_evap_datapoints,
    n_used_evap_datapoints = n_used_evap_datapoints,
    evap_class_id = evap_dat$evap_class_id,
    n_evap_classes = n_evap_classes,
    evap_time = evap_dat$time,
    measured_mass = evap_dat$measured_mass,
    equilibrium_mass = array(evap_class_dat$equilibrium_mass),
    initial_mass = array(evap_class_dat$initial_mass))

flags_list <- list(
    empirical_salt = empirical_salt,
    mass_prior_check = FALSE,
    shared_E_a = shared_E_a)


full_data_list <- c(
    observation_data_list,
    experiment_data_list,
    titer_data_list,
    evap_data_list,
    flags_list)


init_val <- "random"

###############################
## Compile, fit, and save model
###############################
## pass stan data and hyperparams
stan_data <- c(
    full_data_list,
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
