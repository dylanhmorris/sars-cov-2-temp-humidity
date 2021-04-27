#!/usr/bin/env Rscript

########################################
## filename: figure-halflives-compare.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot half-lives for the measured conditions,
## comparing measured versus modeled
#######################################


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(magick))
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
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
evap_dat_path <- args[2]

hl_names <- c(
    "Direct estimate\n",
    "Main mechanistic model\n",
    "Model with directly-\nmeasured concentration")

hl_paths <- c(
    args[3],
    args[4],
    args[5])

outpath <- args[6]

## read data / style files
cat("reading data and images...\n")
dat <- read_data_for_plotting(data_path)

pla_dat <- read_data_for_plotting(data_path)


hl_plots <- list()

for ( hl_ind in 1:length(hl_names) ) {
    hl_type <- hl_names[hl_ind]
    hl_path <- hl_paths[hl_ind]

    chains <- readRDS(hl_path)

##################################################
## calculate posterior draws for half lives
##################################################
    ## get needed draws and add human readable names

    cat("reading chains...\n")
    hl_draws <- chains %>%
        spread_draws(decay_rate[experiment_id]) %>%
        add_titer_metadata(dat, "experiment_id") %>%
        mutate(half_life = to_half_life(decay_rate))

###################################
## plot half-lives as eyeplots
###################################

    temps <- unique(hl_draws$temperature) %>% sort(decreasing = TRUE)
    hums <- unique(hl_draws$humidity) %>% sort()

    halflives <- hl_draws %>%
        ggplot(aes(x = humidity,
                   y = half_life,
                   fill = temperature))

    for(temp in temps){
        halflives <- halflives +
            stat_summary(
                aes(color = temperature,
                    group = temperature),
                fun = "median",
                geom = "line",
                size = 2.5,
                data = hl_draws %>% filter(temperature == temp)) +
            stat_eye(color = "black",
                     alpha = 0.65,
                     .width = 0.95,
                     shape = 21,
                     interval_alpha = 0,
                     position = "identity",
                     data = hl_draws %>% filter(temperature == temp)) +
            stat_slab(
                aes(color = temperature),
                alpha = 1,
                side = "both",
                size = 1,
                fill = NA,
                position = "identity",
                data = hl_draws %>% filter(temperature == temp))
    }

    for(temp in temps){
        halflives <- halflives +
            stat_summary(
                aes(x = humidity,
                    y = half_life,
                    fill = temperature),
                fun = "median",
                color = "black",
                geom = "point",
                shape = 21,
                stroke = 1.25,
                size = 6,
                data = hl_draws %>%
                    filter(temperature == temp))
    }
    
    halflives <- halflives +
        scale_color_temperature(
            guide = "none") +
        scale_fill_temperature(breaks = temps,
                               labels = paste0(temps, "\u00B0C"),
                               guide = guide_legend()) +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(breaks = seq(0, 120, 20),
                           labels = c(seq(0, 100, 20), "evaporation phase")) +
        coord_cartesian(ylim = c(5/10, 200),
                        expand = 0) +
        labs(subtitle = hl_type)

    ## styling: no facet labels because is background plot

    halflives <- halflives +
        theme_project() +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.placement = "outside",
            strip.switch.pad.grid = unit("0.5", "in"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "none")

    if (hl_type == "Direct estimate\n") {
        halflives <- halflives +
            ylab("half-life (hours)") +
            theme(
                legend.key.size = unit(2, "lines"),
                legend.position = c(0.8, 0.85),
                legend.margin = margin(t = -10,
                                       l = 10,
                                       r = 10,
                                       b = 10),
                legend.title = element_blank(),
                legend.box.background = element_rect(color = "black",
                                                     size = 2))
    } else {
        halflives <- halflives + ylab("")
    }

    if (hl_type == "Main mechanistic model\n") {
        halflives <- halflives + 
            xlab("relative humidity (%)")
    } else {
        halflives <- halflives +
            xlab("")
    }
    
    hl_plots[hl_type] <- list(halflives)
}

full_plot <- plot_grid(
    hl_plots[[ hl_names[1] ]],    
    hl_plots[[ hl_names[2] ]],
    hl_plots[[ hl_names[3] ]],
    ncol = 3,
    align = "tblr")

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')

save_plot(
    outpath = outpath,
    fig = full_plot,
    base_height = 10,
    base_asp = 2,
    limitsize = FALSE)

warnings()
