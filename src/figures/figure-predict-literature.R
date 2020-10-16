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

pi_dodge <- position_dodge(width = 0.25)
interval_alpha <- 1


temps_medium <- tibble(temperature = seq(0, 40, length.out = 20))
hums_long = tibble(humidity = seq(0.001, 1, length.out = 100))


scale_stroke <- function(...) {

    return ( scale_discrete_manual(
    aesthetics = "stroke",
    values = unlist(list(
        "RML" = 1.5,
        "literature" = 1)),
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
    mutate(source = "literature")

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
    add_titer_metadata(read_csv(sars_mers_dat_path),
                       "experiment_id") %>%
    mutate(source = "RML")

sars_mers_transient_draws <- sars_mers_hl_chains %>%
    spread_draws(transient_log_half_life[transient_class_id]) %>%
    add_titer_metadata(sars_mers_dat,
                       "transient_class_id") %>%
    mutate(source = "RML")

hl_draws <- lit_hl_draws %>%
    bind_rows(pla_hl_draws) %>%
    bind_rows(sars_mers_hl_draws) %>%
    arrange((virus))

transient_hl_draws <- pla_transient_draws %>%
    bind_rows(sars_mers_transient_draws) %>%
    arrange((virus))

chosen_draws <- sample(1:max(pla_mech_draws['.draw']), n_lines)

##########################
## predict half-lives
## as a function of T/RH
##########################

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


medians <- hl_draws %>% group_by(experiment_id,
                                 temperature,
                                 humidity,
                                 source) %>%
    summarise(log_half_life = median(log_half_life),
              virus = virus[1])

literature_medians <- medians %>%
    filter(source == "literature") %>%
    arrange(desc(virus))

RML_medians <- medians %>%
    filter(source == "RML") %>%
    arrange(desc(virus))

transient_medians <- transient_hl_draws %>%
    group_by(experiment_id,
             temperature,
             humidity,
             source) %>%
    summarise(transient_log_half_life =
                  median(transient_log_half_life),
              virus = factor(virus[1],
                             ordered = TRUE,
                             levels = levels(virus))) %>%
    arrange(desc(virus))

print(transient_medians)


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
        aes(group = interaction(experiment_id, virus, source),
            y = exp(log_half_life),
            x = humidity,
            shape = virus),
        data = hl_draws %>%
            filter(temperature <= 40 &
                   !is.na(humidity) &
                   !is.na(temperature)) %>%
            arrange(virus),
        color = "black",
        position = pi_dodge,
       point_alpha = 0,
        interval_alpha = interval_alpha) +
    stat_pointinterval(
        aes(group = interaction(experiment_id, virus, source),
            y = exp(transient_log_half_life),
            x = 100,
            fill = temperature,
            shape = virus),
        data = transient_hl_draws,
        color = "black",
        position = pi_dodge,
        interval_alpha = interval_alpha,
        point_alpha = 0) +
    geom_point(
        aes(group = interaction(experiment_id, virus, source),
            y = exp(log_half_life),
            x = humidity,
            shape = virus),
        data = literature_medians,
        color = "black",
        position = pi_dodge,
        size = literature_size) +
    geom_point(
        aes(group = interaction(experiment_id, virus, source),
            y = exp(log_half_life),
            x = humidity,
            shape = virus),
        data = RML_medians,
        color = "black",
        position = pi_dodge,
        size = RML_size) +
    geom_point(
        aes(group = interaction(experiment_id, virus, source),
            y = exp(transient_log_half_life),
            x = 100,
            fill = temperature,
            shape = virus),
        data = transient_medians,
        color = "black",
        position = pi_dodge,
        size = RML_size) +
    scale_fill_temperature() +
    scale_color_temperature(guide = FALSE) +
    scale_shape_virus(
        limits = c("SARS-CoV-2", "SARS-CoV-1", "MERS-CoV"),
        guide = guide_legend(reverse = TRUE,
                                           override.aes = list(
                                               fill = "black",
                                               size = c(
                                                   literature_size,
                                                   literature_size,
                                                   literature_size)))) +
    scale_y_continuous(trans = "log10") +
    xlab("relative humidity (%)") +
    ylab("half-life (hours)") +
    coord_cartesian(ylim = c(1e-1, 1e2),
                    xlim = c(0, 102))

save_plot(outpath = out_path,
          fig = lit_panel,
          base_height = 8,
          base_asp = 2)

