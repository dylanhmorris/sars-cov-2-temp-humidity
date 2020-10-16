#!/usr/bin/env Rscript

########################################
## filename: figure-sars-mers-regression.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot fit of simple regression to SARS-CoV-1
## and MERS-CoV data
#######################################


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(virusenv))


#' predict_log_titers
#'
#' calculates predicted log titer
#' at a given time, given intercept,
#' transient (evaporation) phase decay rate,
#' quasi-equilibrium decay rate, and drying time
#' 
#' @param time time of titer observation
#' @param intercept titer at time = 0 in log10 TCID50/0.1mL
#' @param transient_decay_rate inactivation rate during the
#' evaporation ("transient") phase
#' @param decay_rate inactivation rate during the quasi-equilibrium
#' ("dry") phase
#' @param drying_time time to reach quasi-equilibrium phase
#' @return predicted titer in log10 TCID50/mL
#' added to a ggplot plot
predict_log_titers <- function(data){
    return( data %>%
            mutate(predicted_log_titer =
                       1 + predict_titers_implicit_evaporation(
                               time,
                               intercept,
                               transient_decay_rate,
                               decay_rate,
                               drying_time)))
}

#################################
## styling
#################################
gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")

titer_ylab <- expression(atop("virus titer",
                              "(TCID"[50] * "/mL media)"))
titer_xlab <- "time since deposition (hours)"
LOD_log10 <- 0.5
LOD <- 10^LOD_log10
detection_linesize = 0.75
n_lines_decay <- 10
line_alpha_decay <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]
ylim <- c(10^(-0.5), 10^5)
xlim <- c(0, 100)

times <- tibble(
    time = seq(xlim[1], xlim[2], length.out = 250))


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
evap_dat_path <- args[2]
results_path <- args[3]
titers_path <- args[4]
outpath <- args[5]

## read data / style files
cat("reading data...\n")
dat <- read_data_for_plotting(data_path)

pla_dat <- read_data_for_plotting(data_path)

evap_dat <- read_csv(evap_dat_path,
                     col_types = cols()) %>%
    rename(mass = measured_mass)


chains <- readRDS(results_path)
decay_chains <- readRDS(results_path)
titer_ests_chains <- readRDS(titers_path)

cat("extracting draws for decay rates / intercepts...\n")
int_draws <- decay_chains %>%
    spread_draws(intercept[titer_id]) %>%
    add_titer_metadata(pla_dat)

cat("extracting decay rates...\n")
decay_draws <- decay_chains %>%
    spread_draws(c(transient_decay_rate,
                   decay_rate)[experiment_id],
                 drying_time[evap_class_id])


cat("extracting positive wells...\n")
pos_wells <- pla_dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

tidy_draws <- decay_draws %>%
    inner_join(int_draws,
               by = c("experiment_id",
                      ".draw"))

cat('extracting titer estimates...\n')

titer_ests_draws <- titer_ests_chains %>%
    spread_draws(sampled_titer[titer_id])

## get human readable names and detectability
titer_ests_draws <- titer_ests_draws %>%
    add_titer_metadata(pla_dat) %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 1) %>%
    filter(material == "Plastic")


titer_ests_draws <- titer_ests_draws %>%
    mutate(
        log10_tcid50 = ifelse(
            detectable,
            sampled_titer + 1,
            LOD_log10))


###################################
## plot panel showing fit of
## regression lines to real data
###################################
cat('plotting regression lines...\n')
## draw n_lines random regression lines
chosen_draws <- sample(1:max(tidy_draws$.draw), n_lines_decay)

func_samples <- tidy_draws %>%
    filter(.draw %in% chosen_draws)

## cross product decay_rates with x (time) values
## and calculate y (titer) values
cat('setting up x values...\n')

to_plot <- func_samples %>%
    select(-time) %>%
    crossing(times)

## adding one to convert to per mL from per 0.1 mL
to_plot <- to_plot %>%
    predict_log_titers() %>%
    mutate(
        predicted_titer = 10^predicted_log_titer)

shape_scale <- scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))


regression_panel <- to_plot %>%
    ggplot() +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    geom_line(aes(
        x = time,
        y = predicted_titer,
        color = virus,
        group = interaction(.draw, titer_id)),
        alpha = line_alpha_decay) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(x = time,
                      y = 10^log10_tcid50,
                      shape = detectable,
                      fill = virus,
                      group = titer_id),
        data = titer_ests_draws,
        point_size = 6,
        size = 7,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    scale_y_log10_mathformat() +
    scale_fill_virus() +
    scale_color_virus() +
    shape_scale +
    facet_wrap(vars(virus)) +
    xlab(titer_xlab) +
    ylab(titer_ylab) +
    theme_project() +
    coord_cartesian(xlim = xlim,
                    ylim = ylim) +
    theme(legend.position = "none") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


## save the plot to outpath
cat('saving figure to ', outpath, '...\n')

save_plot(
    outpath = outpath,
    fig = regression_panel,
    base_height = 10,
    base_asp = 1.25,
    limitsize = FALSE)
