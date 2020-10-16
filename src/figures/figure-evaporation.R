#!/usr/bin/env Rscript

########################################
## filename: figure-evaporation.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot measured evaporation and estimated
## evaporation rate
#######################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(virusenv))

##########
# styling
##########
n_lines <- 100
line_alpha = 0.1
fill_blue <- "#00b7db"

#################
## data read-in
#################
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
    args[1] <- '../../dat/cleaned/plastic-data.csv'
    args[2] <- '../../dat/cleaned/dmem-evaporation-data.csv'
    args[3] <- '../../out/mcmc_chains/plastic-chains.Rds'
    args[4] <- '../../out/figures/figure-evaporation.pdf'    
}

pla_dat_path <- args[1]
evap_dat_path <- args[2]
chains_path <- args[3]
out_path <- args[4]

pla_dat <- read_data_for_plotting(pla_dat_path)
evap_dat <- read_csv(evap_dat_path,
                     col_types = cols()) %>%
    rename(mass = measured_mass)
chains <- readRDS(chains_path)

draws <- chains %>% spread_draws(
                        beta[evap_class_id]) %>%
    add_titer_metadata(pla_dat, "evap_class_id") %>%
    inner_join(evap_dat %>% distinct(temperature,
                                     humidity,
                                     .keep_all = TRUE) %>%
               select(temperature,
                      humidity,
                      initial_mass,
                      equilibrium_mass),
               by = c("temperature", "humidity"))

chosen_draws <- sample(1:max(draws$.draw), n_lines)

times <- tibble(
    time = seq(0, 60, length.out = 1000))

plot_draws <- draws %>%
    filter(.draw %in% chosen_draws) %>%
    select(-time) %>%
    crossing(times) %>%
    mutate(
        pred_mass = initial_mass - beta * time,
        mass = ifelse(pred_mass >= equilibrium_mass,
                      pred_mass,
                      equilibrium_mass))

fig <- plot_draws %>%
    ggplot(aes(
        x = time,
        y = mass)) +
    geom_line(aes(group = .draw),
              alpha = line_alpha,
              color = "black",
              size = 1) +
    geom_point(shape = 21,
               fill = fill_blue,
               alpha = 0.7,
               size = 2.5,
               stroke = 0.75,
               data = evap_dat) +
    geom_point(
        aes(y = equilibrium_mass),
        x = 24,
        shape = 22,
        fill = fill_blue,
        alpha = 1,
        size = 5,
        stroke = 0.75,
        data = evap_dat %>%
            distinct(temperature,
                     humidity,
                     .keep_all = TRUE)) +
    facet_grid(vars(humidity),
               vars(temperature)) +
    scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
    coord_cartesian(ylim = c(-5, 55),
                    xlim = c(-2.5, 25),
                    expand = FALSE)
######
## style and save figure
######
labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.05, 0.5))

left_margin <- theme(
    plot.margin = margin(b = 3, t = 1, l = 1, r = 3, unit = "cm"))

fig <- fig +
    theme_project() +
    xlab("time since deposition (hours)") +
    ylab("mass (mg)") +
    labs(tag = "relative humidity (%)",
         subtitle = "temperature (\u00B0C)") +
    labeled_panel_theme +
    left_margin +
    theme(legend.position = "none",
          panel.spacing.y = unit(5, "lines"))


## save the plot to outpath
cat('saving figure to ', out_path, '...\n')
save_plot(out_path,
          fig,
          base_height = 15,
          base_asp = 1.2)

warnings()
