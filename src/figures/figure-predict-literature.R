#!/usr/bin/env Rscript

########################################
## filename: figure-predict-literature.R
## author: Dylan Morris <dhmorris@princeton.edu>
## predict exact literature half-lives
## using SARS-CoV-2 mechanistic parameters
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

theme_set(theme_project())

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
    args[1] <- '../../dat/cleaned/plastic-data.csv'
    args[2] <- '../../dat/cleaned/dmem-evaporation-data.csv'
    args[3] <- '../../dat/cleaned/literature-data.csv'
    args[4] <- '../../dat/cleaned/sars-mers-data.csv'

    args[5] <- '../../out/mcmc_chains/plastic-chains.Rds'
    args[6] <- '../../out/mcmc_chains/literature-chains.Rds'
    args[7] <- '../../out/mcmc_chains/sars-mers-chains.Rds'

    args[8] <- '../../out/mcmc_chains/mechanistic-evaporation-modeled-chains.Rds'

    args[9] <- '../../out/figures/figure-predict-literature.pdf'    
}

#################
## styling setup
#################
model_pred_ylim <- 10^c(-1, 2)
ERH <- 0.45
init_frac <- 0.01

n_lines <- 100

literature_size <- 5
RML_size <- 8
    
fixed_lim <- c(1e-1, 1e2)
stroke <- 2
int_big <- 3
int_small <- 1.5
fatten_point = 3
point_size <- 6

pi_dodge <- position_dodge(width = 0.5)
interval_alpha <- 1

transient_humidity <- 100


temps_medium <- tibble(temperature = seq(0, 40, length.out = 20))
hums_long = tibble(humidity = seq(0.001, 1, length.out = 100))

scale_size <- function(point_size = 2,
                       ...) {
    
    return ( scale_discrete_manual(
    aesthetics = "size",
    values = unlist(list(
        "RML" = 1.25  * point_size,
        "literature" = point_size)),
    guide = FALSE,
    ...))
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



cat("Reading data...\n")
pla_dat <- read_data_for_plotting(pla_dat_path)

evap_dat <- read_csv(evap_dat_path,
                     col_types = cols())

lit_dat <- read_data_for_plotting(lit_dat_path) %>%
    mutate(
        super = (is.na(humidity) | humidity > ERH))

sars_mers_dat <- read_data_for_plotting(sars_mers_dat_path)

## read in MCMC chains
cat("Reading MCMC chains...\n")
pla_hl_chains <- readRDS(pla_hl_path)
lit_hl_chains <- readRDS(lit_hl_path)
sars_mers_hl_chains <- readRDS(sars_mers_hl_path)

pla_mech_chains <- readRDS(pla_results_mech_path)

pla_mech_draws <- pla_mech_chains %>%
    spread_draws(E_a_dry,
                 E_a_wet,
                 log_A_dry,
                 log_A_wet,
                 log_hl_20_dry,
                 log_hl_20_wet,
                 log_shape_parameter_concavity,
                 log_shape_parameter_steepness)

## extract MCMC draws into tidy tibbles
cat("Extracting posterior draws...\n")
lit_hl_draws <- lit_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(lit_dat, "experiment_id") %>%
    mutate(data_source = "literature")

pla_hl_draws <- pla_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(pla_dat,
                       "experiment_id") %>%
    mutate(data_source = "RML")

pla_transient_draws <- pla_hl_chains %>%
    spread_draws(transient_log_half_life[transient_class_id]) %>%
    add_titer_metadata(pla_dat,
              "transient_class_id") %>%
    mutate(data_source = "RML")

sars_mers_hl_draws <- readRDS(sars_mers_hl_path) %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(read_csv(sars_mers_dat_path,
                                col_types = cols()),
                       "experiment_id") %>%
    mutate(data_source = "RML")

sars_mers_transient_draws <- sars_mers_hl_chains %>%
    spread_draws(transient_log_half_life[transient_class_id]) %>%
    add_titer_metadata(sars_mers_dat,
                       "transient_class_id") %>%
    mutate(data_source = "RML")

transient_hl_draws <- pla_transient_draws %>%
    bind_rows(sars_mers_transient_draws) %>%
    mutate(log_half_life = transient_log_half_life,
           humidity = transient_humidity) %>%
    ungroup() 

hl_draws <- lit_hl_draws %>%
    bind_rows(pla_hl_draws) %>%
    bind_rows(sars_mers_hl_draws) %>%
    arrange((virus)) %>%
    group_by(experiment_id,
             temperature,
             humidity,
             virus,
             study_author_1,
             data_source) %>%
    mutate(point_id = cur_group_id()) %>%
    ungroup()


chosen_draws <- sample(1:max(pla_mech_draws['.draw']), n_lines)

##########################
## predict half-lives
## as a function of T/RH
##########################

if (!return_csv) {
    cat("Making mechanistic predictions ",
        "(this may take a while)...\n")

    plot_draws_mech <- pla_mech_draws %>%
        filter(.draw %in% chosen_draws) %>%
        crossing(temps_medium) %>%
        crossing(hums_long) %>%
        rowwise() %>%
        mutate(
            log_conc = predict_log_concentration_factor(
                humidity,
                log_shape_parameter_concavity,
                log_shape_parameter_steepness,
                init_frac),
        
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
}

pred_actual_draws <- hl_draws %>%
    filter(!is.na(temperature) &
           !is.na(humidity)) %>%
    inner_join(pla_mech_draws,
               by = ".draw") %>%
    mutate(
        log_conc = predict_log_concentration_factor(
            humidity / 100,
            log_shape_parameter_concavity,
            log_shape_parameter_steepness,
            init_frac),
        
        true_hl = exp(log_half_life),
        
        pred_hl = mechanistic_decay_rate(
            humidity / 100,
            temperature %>% to_Kelvin(),
            log_conc,
            ERH,
            log_A_dry,
            E_a_dry,
            log_A_wet,
            E_a_wet) %>%
            to_half_life())

meds <- pred_actual_draws %>%
    group_by(point_id,
             virus,
             temperature,
             humidity,
             material,
             data_source,
             study_author_1) %>%
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
            temperature,
            humidity,
            data_source,
            study_author_1)


if( return_csv ){
    cat("Making table...\n")
    write_csv(meds %>% select(-point_id),
              out_path)
    
} else {
    cat("Making plot...\n")
    lit_panel <- plot_draws_mech %>%
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
            color = "grey",
            alpha = 1) +
        stat_pointinterval(
            aes(group = point_id,
                y = exp(log_half_life),
                x = humidity,
                shape = virus),
            data = hl_draws %>% filter(data_source == "RML"),
            color = "black",
            position = pi_dodge,
            point_alpha = 1,
            point_size = RML_size,
            interval_alpha = interval_alpha) +
        stat_pointinterval(
            aes(group = point_id,
                y = exp(log_half_life),
                x = humidity,
                shape = virus),
            data = hl_draws %>% filter(data_source == "literature"),
            color = "black",
            position = pi_dodge,
            point_alpha = 1,
            point_size = literature_size,
            interval_alpha = interval_alpha) +
        scale_fill_temperature() +
        scale_color_temperature(guide = FALSE) +
        scale_shape_virus(
            limits = c("SARS-CoV-2", "SARS-CoV-1", "MERS-CoV"),
            guide = FALSE) +
        scale_y_continuous(trans = "log10") +
        xlab("relative humidity (%)") +
        ylab("half-life (hours)") +
        coord_cartesian(ylim = c(1e-1, 1e2),
                        xlim = c(0, 102))

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
            alpha = 0.8) +
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
                                                point_size))))

    cat("Assembling plot from panels...\n")

    temp_leg <- get_legend(
        lit_panel + theme(legend.box.margin = margin(t = 0,
                                                     b = 0,
                                                     l = -50,
                                                     r = 5)))
    virus_leg <- get_legend(
        pred_actual_panel + theme(legend.box.margin = margin(t = 0,
                                                             b = 0,
                                                             l = -50,
                                                             r = 5)))

    leg <- plot_grid(
        temp_leg,
        virus_leg,
        nrow = 3,
        align = "v")


    full_fig <- plot_grid(lit_panel +
                          theme(legend.position = "none"),
                          pred_actual_panel +
                          theme(legend.position = "none"),
                          leg,
                          ncol = 3,
                          align = "hv",
                          axis = "tblr",
                          label_size = 30,
                          labels = c("a", "b", ""),
                          rel_widths = c(1, 1, 0.5))


    save_plot(outpath = out_path,
              fig = full_fig,
              base_height = 8,
              base_asp = 2.6)

}
