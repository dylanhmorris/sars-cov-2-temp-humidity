#!/usr/bin/env Rscript

########################################
## filename: figure-extrapolation.R
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
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(svglite))

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

    args[9] <- '../../out/figures/figure-extrapolation.pdf'    
}


#################
## styling setup
#################
model_pred_ylim <- c(1e-1, 5e2)
ERH <- 0.45
n_lines <- 100

point_size <- 10

    
fixed_lim <- c(1e-4, 1e4)
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
## calculation functions for each panel
###############

#' calc_relative
#'
#' calculates relative prediction
#'
#' @param half_life_chains chains containing estimated half-lives
#' @param mechanistic_draws tidy draws of parameter estimates
#' for mechanistic model
#' @param data_table data corresponding to half-live estimates
calc_relative <- function(half_life_chains, mechanistic_draws, data_table){
    cat("Calculating within-experiment literature predictions...\n")

    lit_ref <- data_table %>%
        distinct(experiment_id, super, .keep_all = TRUE) %>%
        group_by(experiment_class_id, super) %>%
        filter(n_distinct(humidity) > 1 |
               n_distinct(temperature) > 1) %>%
        ungroup() %>%
        mutate(t_dist = abs(temperature - 20),
               h_dist = abs(humidity - ERH * 100)) %>%
        group_by(experiment_class_id, super) %>%
        arrange(t_dist,
                h_dist,
                .by_group = TRUE) %>%
        mutate(ref_exp_id = first(experiment_id))
    
    unique_refs <- lit_ref %>%
        select(ref_exp_id) %>%
        unique()

    expected_points <- dim(lit_ref)[1] - dim(unique_refs)[1]
    cat("Expected number of points: ", expected_points, "\n")


    lit_filt <- data_table %>%
        distinct(experiment_id, .keep_all = TRUE) %>%
        select(ref_temp = temperature,
               ref_hum = humidity,
               ref_exp_id = experiment_id)

    lit_hum <- lit_ref %>%
        inner_join(lit_filt,
                   by = "ref_exp_id") %>%
        ungroup() %>%
        distinct(experiment_id, .keep_all = TRUE)

    experiments_hum <- unique(lit_hum$experiment_id)
    refs_hum <- unique(lit_hum$ref_exp_id)

    hum_lit_draws <- half_life_chains %>%
        spread_draws(log_half_life[experiment_id]) %>%
        inner_join(lit_hum, "experiment_id") %>%
        ungroup()

    hum_lit_draws <- hum_lit_draws %>%
        filter(experiment_id %in% experiments_hum) %>%
        inner_join(hum_lit_draws %>%
                   filter(experiment_id %in% refs_hum) %>%
                   select(ref_exp_id = experiment_id,
                          .draw,
                          ref_log_half_life = log_half_life),
                   by = c("ref_exp_id", ".draw"))


    to_plot <- hum_lit_draws %>%
        inner_join(mechanistic_draws,
                   by = ".draw") %>%
        mutate(
            ref_hl = exp(ref_log_half_life),
            
            true_hl = exp(log_half_life),
                
            E_a_to_use = ifelse(super,
                                E_a_wet,
                                E_a_dry),

            humidity_to_use = replace_na(
                humidity, ERH * 100 + 1),  ## assume equal humidity
            ref_hum_to_use = replace_na(
                ref_hum, ERH * 100 + 1),    ## if not given
        
            pred_hl = half_life_hum(
                ref_hl,
                temperature %>% to_Kelvin(),
                ref_temp %>% to_Kelvin(),
                humidity_to_use / 100,
                ref_hum_to_use / 100,
                ERH,
                E_a_dry,
                E_a_wet,
                log_shape_parameter_concavity)) %>%
    ungroup()

    meds <- to_plot %>%
        filter(experiment_id != ref_exp_id) %>%
        group_by(experiment_id,
                 virus,
                 temperature,
                 humidity,
                 ref_temp,
                 ref_hum,
                 study_author_1,
                 study_doi) %>%
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
            q84_pred = quantile(pred_hl, 0.84)) %>%
        ungroup() %>%
        arrange(desc(virus),
                study_author_1,
                temperature,
                humidity)

    actual_points <- dim(meds)[1]
    cat("Actual number of points", actual_points, "\n")

    if(! (expected_points == actual_points) )
        stop("expected and actual points for predicted/actual plot differ")

    return( meds )
}

#' calc_conc_points
#'
#' calculate measured half-lives from direct regression
#' estimates and where they should fall given
#' concentration factors
calc_conc_points <- function(half_life_draws,
                             mechanistic_draws,
                             evaporation_data,
                             ERH) {
    conc_points <- half_life_draws %>%
        inner_join(mechanistic_draws %>%
                   select(.draw,
                          initial_mass_fraction_solute),
                   by = ".draw") %>%
        inner_join(evaporation_data %>% distinct(temperature,
                                                 humidity,
                                                 .keep_all = TRUE) %>%
                   select(temperature,
                          humidity,
                          initial_mass,
                          equilibrium_mass),
                   by = c("temperature",
                          "humidity")) %>%
        mutate(
            log_concentration_factor = mass_change_to_log_concentration_factor(
                equilibrium_mass,
                initial_mass,
                initial_mass_fraction_solute)) %>%
        filter(humidity > ERH)

    conc_points <- conc_points %>%
        group_by(experiment_id) %>%
        summarise(
            med_conc = exp(median(log_concentration_factor)),
            med_hl = exp(median(log_half_life))) %>%
        inner_join(conc_points, by = "experiment_id")

    return (conc_points)
}


#' calc_conc_to_plot
#'
#' predict halflives as a function
#' of concentration factor from the
#' mechanistic model
calc_conc_to_plot <- function(mechanistic_draws,
                              chosen_draw_ids,
                              temperatures_to_plot,
                              concentrations_to_plot) {
    draws_to_plot <- mechanistic_draws %>%
        filter(.draw %in% chosen_draw_ids) %>%
        crossing(temperatures_to_plot) %>%
        crossing(concentrations_to_plot) %>%
        rowwise() %>%
        mutate(
            log_conc = log(concentration_factor),
            
            decay_rate = mechanistic_decay_rate(
                ERH + 0.05,
                temperature %>% to_Kelvin(),
                log_conc,
                ERH,
                log_A_dry,
                E_a_dry,
                log_A_wet,
                E_a_wet),
            
            half_life_model = decay_rate %>% to_half_life()
        ) %>%
        group_by(.draw, temperature) %>%
        mutate(line_id = cur_group_id()) %>%
        ungroup()

    return( draws_to_plot )
}

#' calc_hum_to_plot
#'
#' predict halflives as a function
#' of relative humidity from the
#' mechanistic model
calc_hum_to_plot <- function(mechanistic_draws,
                             chosen_draw_ids,
                             temperatures_to_plot,
                             humidities_to_plot) {

    draws_to_plot <- mechanistic_draws %>%
        filter(.draw %in% chosen_draw_ids) %>%
        crossing(temperatures_to_plot) %>%
        crossing(humidities_to_plot) %>%
        rowwise() %>%
        mutate(
            log_conc = predict_log_concentration_factor(
                humidity,
                log_shape_parameter_concavity,
                log_shape_parameter_steepness,
                initial_mass_fraction_solute),
        
            decay_rate = mechanistic_decay_rate(
                humidity,
                temperature %>% to_Kelvin(),
                log_conc,
                ERH,
                log_A_dry,
                E_a_dry,
                log_A_wet,
                E_a_wet),
        
            half_life_model = decay_rate %>% to_half_life()
        ) %>%
        group_by(.draw, temperature) %>%
        mutate(line_id = cur_group_id()) %>%
        ungroup()

    return( draws_to_plot )
}


#' matrix_predictions
#' predict half-lives across a matrix
#' of temperatures and relative humidities
matrix_predictions <- function(mechanistic_draws) {
    cat("Making matrix predictions (this may take a while)...\n")
    mat <- tibble(
        temperature = seq(0, 40, length.out = 100)) %>%
        crossing(tibble(
            humidity = seq(0, 100, length.out = 100))) %>%
        crossing(mechanistic_draws) %>%
        mutate(
            log_conc = predict_log_concentration_factor(
                humidity / 100,
                log_shape_parameter_concavity,
                log_shape_parameter_steepness,
                initial_mass_fraction_solute),

            half_life = mechanistic_decay_rate(
                humidity / 100,
                temperature %>% to_Kelvin(),
                log_conc,
                ERH,
                log_A_dry,
                E_a_dry,
                log_A_wet,
                E_a_wet) %>%
                to_half_life()) %>%
    group_by(temperature, humidity) %>%
    summarise(median_half_life = median(half_life))

    return( mat )
}

#' prep_halflives_for_heatmap
#'
#' prepare halflives from the literature
#' and from own experiments for plotting
#' on top of T/RH heatamp
prep_halflives_for_heatmap <- function(merged_draws,
                                       own_main_data,
                                       own_sars_mers_data,
                                       literature_data) {
    hl_to_plot <- merged_draws %>%
        group_by(temperature, humidity, virus, experiment_id, material, source) %>%
        summarise(median_half_life = median(exp(log_half_life)))


    hl_to_plot <- hl_to_plot %>%
        ungroup() %>%
        filter(!is.na(temperature) &
               !is.na(humidity))

    hl_to_plot <- hl_to_plot %>%
        group_by(temperature, humidity) %>%
        mutate(rn = row_number()) %>%
        ungroup()

    hl_to_plot <- hl_to_plot %>%      
        mutate(temperature_nudged = cycle_nudge_y(temperature,
                                                  rn,
                                                  2),
               humidity_nudged = cycle_nudge_x(humidity,
                                               rn,
                                               5)) %>%
        ungroup()

    expected_points_lit <- dim(
        literature_data %>% filter(!is.na(temperature) &
                                   !is.na(humidity)) %>%
        distinct(temperature,
                 humidity,
                 virus,
                 material,
                 study_author_1))[1]
    
    expected_points_main <- dim(
        own_main_data %>% filter(!is.na(temperature) &
                                   !is.na(humidity)) %>%
        distinct(temperature,
                 humidity,
                 virus,
                 material))[1]

    expected_points_sars_mers <- dim(
        own_sars_mers_data %>% filter(!is.na(temperature) &
                                      !is.na(humidity)) %>%
        distinct(temperature,
                 humidity,
                 virus,
                 material))[1]

    expected_points <- (expected_points_lit +
                        expected_points_sars_mers +
                        expected_points_main)
    
    cat("expected number of points:", expected_points, "\n")

    actual_points <- dim(hl_to_plot)[1]
    cat("Actual number of points", actual_points, "\n")
                                     
    
    if(actual_points != expected_points)
        stop("expected and actual points for heatmap differ")
    
    return( hl_to_plot )
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

pla_results_mech_path <- args[8]

out_path <- args[9]

return_csv <- grepl(".csv", out_path)
if(return_csv){
    cat("Making csv data file", out_path, "...\n")
} else {
    cat("Making figure ", out_path, "...\n")
}


## read in data files
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
pla_mech_chains <- readRDS(pla_results_mech_path)

## extract MCMC draws into tidy tibbles
cat("Extracting posterior draws...\n")
pla_hl_draws <- pla_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(pla_dat,
                       "experiment_id") %>%
    mutate(source = "RML")

pla_transient_draws <- pla_hl_chains %>%
    spread_draws(transient_log_half_life[transient_class_id]) %>%
    add_titer_metadata(pla_dat,
              "transient_class_id") %>%
    mutate(source = "RML")

sars_mers_hl_draws <- readRDS(sars_mers_hl_path) %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(read_csv(sars_mers_dat_path,
                                col_types = cols()),
                       "experiment_id") %>%
    mutate(source = "RML")

sars_mers_transient_draws <- sars_mers_hl_chains %>%
    spread_draws(transient_log_half_life[transient_class_id]) %>%
    add_titer_metadata(sars_mers_dat,
                       "transient_class_id") %>%
    mutate(source = "RML")

lit_hl_draws <- lit_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    ungroup() %>%
    add_titer_metadata(lit_dat, "experiment_id") %>%
    mutate(source = "literature")

hl_draws <- lit_hl_draws %>%
    bind_rows(pla_hl_draws) %>%
    bind_rows(sars_mers_hl_draws)

transient_hl_draws <- pla_transient_draws %>%
    bind_rows(sars_mers_transient_draws)

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

##########################
## do calculations
##########################
## within-study prediction (normalized) calculations
meds <- calc_relative(lit_hl_chains,
                      pla_mech_draws,
                      lit_dat)

if(!return_csv){
    ## concentration factor calculations
    cat("Calculating predictions for half-life as a",
        "function of concentration factor...\n")

    plot_draws_conc <- calc_conc_to_plot(
        pla_mech_draws,
        chosen_draws,
        temps_short,
        concs_long)

    conc_points <- calc_conc_points(
        hl_draws,
        pla_mech_draws,
        evap_dat,
        ERH * 100)


    ## humidity calculations
    cat("Calculating predictions for half-life as a",
        "function of relative humidity...\n")

    plot_draws_hum <- calc_hum_to_plot(
        pla_mech_draws,
        chosen_draws,
        temps_short,
        hums_long)
    

    ## heatmap calculations
    cat("Calculating heatmap with literature data...\n")

    mat <- matrix_predictions(pla_mech_draws)

    cat("preparing literature half-lives for heatmap superimposition...\n")
    hl_to_plot <- prep_halflives_for_heatmap(hl_draws,
                                             pla_dat,
                                             sars_mers_dat,
                                             lit_dat)
}

######################
## do plotting / table manufacture
######################
make_plot <- function(){                      
    cat("plotting modeled humidity/concentration panel...\n")
    hum_panel <- plot_draws_hum %>%
        ggplot(aes(
            x = humidity * 100,
            y = half_life_model,
            color = temperature,
            fill = temperature,
            group = line_id)) +
        geom_line(
            alpha = 0.1) +
        geom_vline(
            xintercept = ERH * 100,
            size = 2,
            linetype = "dashed") +
        stat_pointinterval(
            aes(group = interaction(temperature, humidity),
                y = exp(log_half_life),
                x = humidity,
                shape = virus),
            interval_size_range = c(int_small, int_big),
            fatten_point = fatten_point,
            data = pla_hl_draws,
            color = "black",
            stroke = stroke) +
        stat_pointinterval(
            aes(group = interaction(temperature, humidity),
                y = exp(transient_log_half_life),
                x = 100,
                shape = virus),
            interval_size_range = c(int_small, int_big),
            fatten_point = fatten_point,
            data = transient_hl_draws %>%
                filter(virus == "SARS-CoV-2"),
            position = position_dodge(width = -5),
            color = "black",
            stroke = stroke) +
        scale_fill_temperature(name = temp_label) +
        scale_color_temperature(guide = FALSE) +
        scale_shape_virus(guide = FALSE) + 
        scale_y_log10_mathformat() +
        xlab(rh_label) +
        ylab(hl_label) +
        coord_cartesian(ylim = model_pred_ylim,
                        xlim = c(0, 105),
                        expand = 0)

    cat("plotting concentration factor panel...\n")
    
    conc_panel <- plot_draws_conc %>%
        ggplot(aes(
            x = concentration_factor,
            y = half_life_model,
            color = temperature,
            fill = temperature,
            group = line_id)) +
        geom_line(
            alpha = 0.1) +
        stat_pointinterval(
            aes(group = experiment_id,
                y = exp(log_half_life),
                x = med_conc),
            data = conc_points,
            color = "black",
            interval_size_range = c(int_small, int_big),
            fatten_point = fatten_point,
            stroke = stroke,
            point_alpha = 1,
            shape = 21) +
        stat_pointinterval(
            aes(group = experiment_id,
                y = med_hl,
                x = exp(log_concentration_factor)),
            data = conc_points,
            color = "black",
            interval_size_range = c(int_small, int_big),
            fatten_point = fatten_point,
            stroke = stroke,
            shape = 21) +
        scale_shape_virus(guide = FALSE) +
        scale_fill_temperature() +
        scale_color_temperature() +
        scale_y_continuous(trans = "log10",
                           breaks = c(0.1, 1, 10, 100)) +
        coord_cartesian(xlim = c(1, 100),
                        ylim = model_pred_ylim,
                        expand = 0) +
        xlab("concentration factor") +
        ylab(hl_label) +
        guides(color = guide_legend(override.aes = list(alpha = 1)))


    cat("plotting heatmap panel...\n")
    heatmap_panel <- mat %>%
        ggplot(aes(
            x = humidity,
            y = temperature,
            fill = median_half_life)) +
        geom_raster(data = mat,
                    interpolate = FALSE,
                    hjust = 0,
                    vjust = 0) +
        geom_point(
            data = hl_to_plot %>%
                arrange(source),
            mapping = aes(
                x = humidity_nudged,
                y = temperature_nudged,                
                shape = virus,
                stroke = source,
                size = source,
                group = experiment_id)) +
        scale_size() +
        scale_stroke() + 
        scale_shape_virus(guide = FALSE) +
        scale_fill_viridis(trans = "log10",
                           direction = -1,
                           name = paste0("median ", hl_label)) +
        coord_cartesian(xlim = c(0, 105),
                        ylim = c(0, 40),
                        expand = 0) +
        ylab(temp_label) +
        xlab("relative humidity (%)")


    cat("Plotting all within-experiment prediction panel...\n")
    pred_actual_panel <- meds %>%
        ggplot(aes(
            x = med_pred,
            y = med_true,
            fill = virus,
            shape = virus,
            group = experiment_id)) +
        geom_abline(intercept = 0,
                    slope = 1,
                    size = stroke) +
        geom_linerange(
            size = int_small * 0.5,
            aes(xmin = q025_pred,
                xmax = q975_pred)) +    
        geom_linerange(
            size = int_big * 0.5,
            aes(xmin = q16_pred,
                xmax = q84_pred)) +
        geom_linerange(
            size = int_small * 0.5,
            aes(ymin = q025_true,
                ymax = q975_true)) +
        geom_linerange(
            size = int_big * 0.5,
            aes(ymin = q16_true,
                ymax = q84_true)) +
        geom_point(
            aes(x = med_pred,
                y = med_true,
                size = virus,),
            stroke = stroke,
            alpha = 0.9) +
        scale_shape_virus() +
        scale_size_virus(base_size = point_size,
                         guide = FALSE) +
        scale_fill_virus() +
        scale_x_log10_mathformat() +
        scale_y_log10_mathformat() +
        coord_cartesian(xlim = fixed_lim,
                        ylim = fixed_lim,
                        expand = 0) +
        ylab("measured half-life (hours)") +
        xlab("predicted half-life (hours)") +
        guides(fill = guide_legend(override.aes = list(
                                       size = c(point_size,
                                                point_size,
                                                point_size,
                                                0.75 * point_size,
                                                0.75 * point_size))))


    ## assemble figure from panels
    temp_leg <- get_legend(
        hum_panel + theme(legend.box.margin = margin(0, 0, 0, 30)))
    virus_leg <- get_legend(
        pred_actual_panel + theme(legend.box.margin = margin(0, 0, 0, 30)))
    hl_leg <- get_legend(
        heatmap_panel + theme(legend.box.margin = margin(0, 0, 0, 30)))

    leg <- plot_grid(
        temp_leg,
        hl_leg,
        virus_leg,
        ncol = 3,
        align = "h")

    aligned_plots <- align_plots(
        hum_panel + theme(legend.position = "none",
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank()),
        conc_panel + theme(legend.position = "none",
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank()),
        heatmap_panel + theme(legend.position = "none"),
        pred_actual_panel + theme(legend.position = "none"),
        align = "hv",
        axis = "tblr",
        greedy = FALSE)

    full_fig <- plot_grid(
        aligned_plots[[1]],
        aligned_plots[[2]],
        aligned_plots[[3]],
        aligned_plots[[4]],
        label_size = 30,
        labels = "auto",
        ncol = 2,
        nrow = 2)

    fig_with_legend <- plot_grid(full_fig,
                                 leg,
                                 ncol = 1,
                                 rel_heights = c(5, 1))

    save_plot(outpath = out_path,
              fig = fig_with_legend,
              base_height = 20,
              base_asp = 5 / 6)
}

make_table <- function(table_to_write,
                       out_path){
    table_to_write %>% write_csv(out_path)
}


## make plot or table
if( return_csv) {
    cat("Making table!\n")
    make_table(meds, out_path)
} else {
    cat("Making plot!\n")
    make_plot()
}
