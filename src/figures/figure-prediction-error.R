#!/usr/bin/env Rscript

########################################
## filename: figure-prediction-error.R
## author: Dylan Morris <dhmorris@princeton.edu>
## extrapolate from mech model to
## unobserved conditions
#######################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(virusenv))

args <- commandArgs(trailingOnly = TRUE)

theme_set(theme_project() +
          theme(panel.border = element_blank()))

if(length(args) < 1) {
    args[1] <- '../../dat/cleaned/plastic-data.csv'
    args[2] <- '../../dat/cleaned/dmem-evaporation-data.csv'

    args[3] <- '../../dat/cleaned/literature-data.csv'
    args[4] <- '../../dat/cleaned/sars-mers-data.csv'
    args[5] <- '../../out/mcmc_chains/plastic-chains.Rds'
    args[6] <- '../../out/mcmc_chains/literature-chains.Rds'
    args[7] <- '../../out/mcmc_chains/sars-mers-chains.Rds'

    args[8] <- '../../out/mcmc_chains/mechanistic-evaporation-modeled-chains.Rds'
    args[9] <- '../../out/mcmc_chains/mechanistic-evaporation-modeled-chains.Rds'

    args[10] <- '../../out/figures/figure-prediction-error.pdf'
    theme_set(theme_minimal())
}


#################
## styling setup
#################
model_pred_ylim <- c(1e-1, 5e2)
ERH <- 0.45
n_lines <- 100

point_size <- 10

    
fixed_lim <- c(1e-1, 1e2)
stroke <- 2
int_big <- 3
int_small <- 1.5
fatten_point = 3


temp_label <- "temperature (\u00B0C)"
rh_label <- "relative humidity (%)"
hl_label <- "half-life (hours)"

scale_stroke <- function(...) {

    return ( scale_discrete_manual(
    aesthetics = "stroke",
    values = unlist(list(
        "RML" = 1.5,
        "literature" = 1)),
    guide = FALSE,
    ...))
}

scale_size <- function(...) {

    return ( scale_discrete_manual(
    aesthetics = "size",
    values = unlist(list(
        "RML" = 1.25  * point_size,
        "literature" = point_size)),
    guide = FALSE,
    ...))
}

compass_rotate_x <- function(compass_point){
    x <- (compass_point %% 2) * (compass_point - 2)
    return( x )
}
compass_rotate_y <- function(compass_point){
    y <- -((compass_point - 1) %% 2) * (compass_point - 3)
    return( y )
}

cycle_nudge_x <- function(x, count, nudge_per_cycle = 1) {
    val <- count - 1
    cycle <- ceiling(val / 4)
    compass_point <- (val %% 4) + 1
    return (x + cycle * compass_rotate_x(compass_point))
}

cycle_nudge_y <- function(y, count, nudge_per_cycle = 1) {
    val <- count - 1
    cycle <- ceiling(val / 4)
    compass_point <- (val %% 4) + 1
    return (y + cycle * compass_rotate_y(compass_point))
}


########################
## variable setup                                     #
########################

temps_long <- tibble(temperature = seq(0, 40, length.out = 100))
temps_short <- tibble(temperature = c(0, 10, 22, 27, 40))

hums_short <- tibble(humidity = c(0.4, 0.65, 0.85, 0.98))
hums_long = tibble(humidity = seq(0, 0.999, length.out = 100))

concs_long = tibble(concentration_factor = seq(1, 100, length.out = 100))

modeled_log_conc <- function(fractional_humidity,
                             initial_fraction,
                             ERH){

    ## avoid throwing error from Tang for
    ## low humidities that we don't
    ## need anyway
    hums_calc <- ifelse(
        fractional_humidity > ERH,
        fractional_humidity,
        0.5)

    tang_mass_fracs <- sapply(hums_calc, tang_mass_frac)
    
    concs <- fraction_change_to_log_concentration_factor(
        tang_mass_fracs,
        initial_fraction)

    return (ifelse(concs > 0,
                   concs,
                   0))
}

#################
## data read-in
#################
pla_dat_path <- args[1]
evap_dat_path <- args[2]
lit_dat_path <- args[3]
sars_mers_dat_path <- args[4]

pla_hl_path <- args[5]
lit_hl_path <- args[6] 
sars_mers_hl_path <- args[7]

pla_results_emp_path <- args[8]
pla_results_mech_path <- args[9]

out_path <- args[10]


cat("Reading data...\n")
pla_dat <- read_data_for_plotting(pla_dat_path)
evap_dat <- read_csv(evap_dat_path,
                     col_types = cols())

lit_dat <- read_data_for_plotting(lit_dat_path) %>%
    mutate(super = (is.na(humidity) | humidity / 100 > ERH))

                       

sars_mers_dat <- read_data_for_plotting(sars_mers_dat_path)

## read in MCMC chains
cat("Reading MCMC chains...\n")
pla_hl_chains <- readRDS(pla_hl_path)
lit_hl_chains <- readRDS(lit_hl_path)
sars_mers_hl_chains <- readRDS(sars_mers_hl_path)

pla_emp_chains <- readRDS(pla_results_emp_path)
pla_mech_chains <- readRDS(pla_results_mech_path)


## extract MCMC draws into tidy tibbles
cat("Extracting posterior draws...\n")
pla_hl_draws <- pla_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(pla_dat,
                       "experiment_id") %>%
    mutate(source = "RML")

sars_mers_hl_draws <- readRDS(sars_mers_hl_path) %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(read_csv(sars_mers_dat_path),
                       "experiment_id") %>%
    mutate(source = "RML")

lit_hl_draws <- lit_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    ungroup() %>%
    add_titer_metadata(lit_dat, "experiment_id") %>%
    mutate(source = "literature")

hl_draws <- lit_hl_draws %>%
    bind_rows(sars_mers_hl_draws)

pla_mech_draws <- pla_mech_chains %>%
    spread_draws(E_a_dry,
                 E_a_wet,
                 log_A_dry,
                 log_A_wet,
                 log_hl_20_dry,
                 log_hl_20_wet,
                 log_shape_parameter_concavity,
                 log_shape_parameter_steepness,
                 initial_mass_fraction_solute,
                 sd_logconc)

chosen_draws <- sample(1:max(pla_mech_draws['.draw']), n_lines)


##############################
## within-study predictions
##############################

cat("Calculating absolute literature predictions...\n")


to_plot <- hl_draws %>%
    filter(!is.na(humidity) &
           !is.na(temperature)) %>%
    inner_join(pla_mech_draws,
               by = ".draw") %>%
    mutate(
        true_hl = exp(log_half_life),
        
        log_conc = predict_log_concentration_factor(
            humidity / 100,
            log_shape_parameter_concavity,
            log_shape_parameter_steepness,
            initial_mass_fraction_solute),

        pred_hl = mechanistic_decay_rate(
            humidity / 100,
            temperature %>% to_Kelvin(),
            log_conc,
            ERH,
            log_A_dry,
            E_a_dry,
            log_A_wet,
            E_a_wet) %>%
            to_half_life()) %>%
    ungroup()


meds <- to_plot %>%
    group_by(virus, temperature, humidity,
             source, study_author_1,
             experiment_id) %>%
    summarise(
        med_true = median(true_hl),
        med_pred = median(pred_hl),
        q025_true = quantile(true_hl, 0.025),
        q975_true = quantile(true_hl, 0.975),
        q025_pred = quantile(pred_hl, 0.025),
        q975_pred = quantile(pred_hl, 0.975),
        q16_true = quantile(true_hl, 0.16),
        q84_true = quantile(true_hl, 0.84),
        q16_pred = quantile(pred_hl, 0.16),
        q84_pred = quantile(pred_hl, 0.84),
        point_id = cur_group_id()) %>%
    ungroup()

######################
## do plotting
######################

print(meds)

cat("Plotting all within-experiment prediction panel...\n")
pred_actual_panel <- meds %>%
    ggplot(aes(
        x = med_pred,
        y = med_true,
        fill = virus,
        shape = virus,
        group = point_id)) +
    geom_abline(intercept = 0,
                slope = 1,
                size = stroke) +
    geom_linerange(
        size = int_small * 0.8,
        aes(xmin = q025_pred,
            xmax = q975_pred)) +    
    geom_linerange(
        size = int_big * 0.8,
        aes(xmin = q16_pred,
            xmax = q84_pred)) +
    geom_linerange(
        size = int_small * 0.8,
        aes(ymin = q025_true,
            ymax = q975_true)) +
    geom_linerange(
        size = int_big * 0.8,
        aes(ymin = q16_true,
            ymax = q84_true)) +
    geom_point(
        aes(x = med_pred,
            y = med_true,
            size = virus,
            shape = virus),
        stroke = stroke,
        alpha = 0.9) +
    scale_shape_virus() +
    scale_size_virus(base_size = point_size,
                     guide = FALSE) +
    scale_fill_virus() +
    scale_x_log10_mathformat() +
    scale_y_log10_mathformat() +
    coord_fixed(xlim = fixed_lim,
                ylim = fixed_lim,
                expand = 0) +
    ylab("measured half-life (hours)") +
    xlab("predicted half-life (hours)") +
    guides(fill = guide_legend(override.aes = list(
                                   size = c(point_size,
                                            point_size,
                                            point_size))))

##############################
## save figure
##############################
cat("Saving to", out_path, "...\n")

save_plot(outpath = out_path,
          fig = pred_actual_panel,
          base_height = 15,
          base_asp = 4/5)

