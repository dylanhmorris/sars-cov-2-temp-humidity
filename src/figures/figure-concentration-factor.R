#!/usr/bin/env Rscript

########################################
## filename: figure-concentration-factor.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## show how we model concentration of salts
## as a function of RH
#######################################

suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(extrafont))

args <- commandArgs(trailingOnly = TRUE)

pla_dat_path <- args[1] 
evap_dat_path <- args[2] 
pla_results_mech_path <- args[3]
out_path <- args[4]

pla_dat <- read_data_for_plotting(pla_dat_path)
evap_dat <- read_csv(evap_dat_path,
                     col_types = cols())

pla_mech_chains <- readRDS(pla_results_mech_path)

ERH <- 0.45

int_big <- 3
int_small <- 1.5
fatten_point = 3

cat("Extracting MCMC draws...\n")

pla_mech_draws <- pla_mech_chains %>%
    spread_draws(E_a_dry,
                 E_a_wet,
                 log_A_dry,
                 log_A_wet,
                 log_hl_20_dry,
                 log_hl_20_wet,
                 log_shape_parameter_concavity,
                 log_shape_parameter_steepness,
                 initial_mass_fraction_solute) %>%
    mutate(model = "plastic measured")


cat("Calculating empirical measurements...\n")
masses <- pla_dat %>%
    distinct(experiment_id, .keep_all = TRUE) %>%
    inner_join(evap_dat %>% distinct(temperature,
                                     humidity,
                                     .keep_all = TRUE) %>%
               select(initial_mass,
                      equilibrium_mass,
                      temperature,
                      humidity),
               by = c("temperature", "humidity")) %>%
    crossing(pla_mech_draws) %>%
    mutate(log_empirical_conc = mass_change_to_log_concentration_factor(
               equilibrium_mass,
               initial_mass,
               initial_mass_fraction_solute),
           
           empirical_concentration_factor = exp(
               log_empirical_conc))
                  

conc_preds <- pla_mech_draws %>%
    crossing(tibble(humidity = seq(40, 100, length.out = 100))) %>%
    mutate(super = (humidity / 100) > ERH)
        

cat("Predicting modeled concentration factors...\n")
conc_pred_draws <- sample(1:max(conc_preds$.draw), 100)


## cat("Tang polynomial prediction...\n")
## tangs <- tibble(humidity = seq(40, 99.9, length.out = 100),
##                 super = humidity / 100 > ERH,
##                 log_conc_tang = ifelse(
##                     super,
##                     fraction_change_to_log_concentration_factor(
##                         sapply(humidity / 100,
##                                tang_mass_frac),
##                         initial_mass_frac),
##                     fraction_change_to_log_concentration_factor(
##                         tang_mass_frac(ERH),
##                         init_mass_frac)),
##                 conc_tang = exp(pmax(0, log_conc_tang)))

cat("Semi-mechanistic function prediction...\n")

print("c_c:")
print(quantile(conc_preds$log_shape_parameter_concavity))
print("\n")
print("c_s:")
print(quantile(conc_preds$log_shape_parameter_steepness))
print("solutes:")
print(quantile(conc_preds$initial_mass_fraction_solute))


conc_preds <- conc_preds %>%
    filter(.draw %in% conc_pred_draws) %>%
    rowwise() %>%
    mutate(
        log_conc = ifelse(
            super,
            predict_log_concentration_factor(
                humidity / 100,
                log_shape_parameter_concavity,
                log_shape_parameter_steepness,
                initial_mass_fraction_solute),
            predict_log_concentration_factor(
                ERH,
                log_shape_parameter_concavity,
                log_shape_parameter_steepness,
                initial_mass_fraction_solute)),
        
        concentration_factor = exp(max(0, log_conc))
    ) %>%
    group_by(.draw) %>%
    mutate(line_id = cur_group_id())


cat("plotting figure...\n")

fig <- masses %>%
    ggplot(aes(
        x = humidity,
        y = empirical_concentration_factor)) +
    geom_vline(xintercept = ERH * 100,
               color = "black",
               linetype = "dashed",
               size = 3) +
    geom_line(
        mapping = aes(
            x = humidity,
            y = concentration_factor,
            group = line_id),
        color = "blue",
        data = conc_preds,
        size = 1,
        alpha = 0.1,
        linetype = "solid") +
    stat_pointinterval(
        aes(fill = temperature,
            group = interaction(temperature, humidity)),
        shape = 21,
        color = "black",
        fatten_point = fatten_point,
        interval_size_range = c(int_small, int_big),
        position = position_dodge(width = 1)) +
    scale_fill_temperature() +
    scale_color_temperature() +
    scale_y_continuous(trans = "log2") +
    ylab("concentration factor") +
    xlab("relative humidity (%)") +
    coord_cartesian(xlim = c(40, 100),
                    ylim = c(1, 128)) +
    theme_project() +
    theme(plot.subtitle = element_text(hjust = 0.5))

cat("Saving plot to", out_path, "...\n")
save_plot(outpath = out_path,
          fig = fig,
          base_height = 10,
          base_asp = 2)

warnings()
