#!/usr/bin/env Rscript

########################################
## filename: figure-halflives-schematic.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot half lives for temp/humid environmental
## conditions analysis alongside schematic
## of the mechanisms of decay
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


#' schematic_curve
#'
#' draws a schematic curve of how
#' half-life should depend on RH
#' 
#' @param y_anchor where on the y axis to anchor the
#' curve, default 0.5
#' @param ... other parameters passed to geom_curve()
#' @return list of geom_curve objects that can be
#' added to a ggplot plot
schematic_curve <- function(y_anchor = 0.5,
                            ...){

    deltas <- c(
        -0.035,
        -0.07,
        -0.28,
        -0.3,
        0.1)
    
    return (
        list(
            geom_curve(
                x = 0,
                y = y_anchor,
                xend = 37,
                yend = y_anchor + deltas[1],
                curvature = 0,
                ...),
            geom_curve(
                x = 37,
                y = y_anchor + deltas[1],
                xend = 40,
                yend = y_anchor + deltas[2],
                curvature = -0.3,
                ...),
            geom_curve(
                x = 40,
                y = y_anchor + deltas[2],
                xend = 47,
                yend = y_anchor + deltas[3],
                curvature = 0.1,
                ...),
            geom_curve(
                x = 47,
                y = y_anchor + deltas[3],
                xend = 50,
                yend =  y_anchor + deltas[4],
                curvature = 0.25,
                ...),
            geom_curve(
                x = 50,
                y =  y_anchor + deltas[4],
                xend = 90,
                yend =  y_anchor + deltas[5],
                angle = 75,
                curvature = 0.25,
                ...)))
}
#' schematic_arrow
#'
#' draws an arrow between parts of the schematic
#' 
#' @param origin_loc vector giving x, y location of the
#' origin for the arrow
#' @param destination_loc vector giving x, y location
#' to which the arrow points
#' @param origin_pad vector giving how far away from the origin to start
#' (in x, y plot units)
#' @param destination_pad vector giving how far away from the destination to
#' stop (in x, y plot units)
#' @param ... other parameters passed to geom_curve()
#' @return customized geom_curve object that can be
#' added to a ggplot plot
schematic_arrow <- function(origin_loc,
                            destination_loc,
                            origin_pad,
                            destination_pad,
                            ...){

    signs = 1 - 2 * (origin_loc > destination_loc)
    
    return ( geom_curve(
        x = origin_loc[1] + signs[1] * origin_pad[1],
        y = origin_loc[2] + signs[2] * origin_pad[2],
        xend = destination_loc[1] - signs[1] * destination_pad[1],
        yend = destination_loc[2] - signs[2] * destination_pad[2],
        lineend = "round",
        arrow = arrow(length = unit(0.03, "npc"),
                      type = "closed"),
        ...) )
}

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
results_path <- args[3]
titers_path <- args[4]
schematic_paths <- args[5:8]
outpath <- args[9]

## read data / style files
cat("reading data and images...\n")
dat <- read_data_for_plotting(data_path)

pla_dat <- read_data_for_plotting(data_path)

evap_dat <- read_csv(evap_dat_path,
                     col_types = cols()) %>%
    rename(mass = measured_mass)


chains <- readRDS(results_path)
decay_chains <- readRDS(results_path)
titer_ests_chains <- readRDS(titers_path)

efflor_img <- image_read(schematic_paths[1])
evap_img <- image_read(schematic_paths[2])
large_img <- image_read(schematic_paths[3])
small_img <- image_read(schematic_paths[4])
theme_family <- "CM Sans"


#################################
## overall plot styling
#################################
set.seed(2342) # reproducible! (since we use random draws)
text_size = 40
chosen_temperature <- 27
chosen_rh <- 85

n_lines_evap <- 10
n_lines_decay <- 10

line_alpha_evap = 0.1
fill_blue <- "#00b7db"
LOD_log10 <- 0.5
LOD <- 10^LOD_log10
detection_linesize = 0.75

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")

titer_ylab <- expression(atop("virus titer",
                              "(TCID"[50] * "/mL media)"))
titer_xlab <- "time since deposition (hours)"

line_alpha_decay <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

ylim_evap <- c(-5, 55)
ylim_titer <- c(10^(-0.5), 1e5)
xlim <- c(-1, 35)
times <- tibble(
    time = seq(0, xlim[2], length.out = 250))


drying_line_size <- 3
drying_line_alpha <- 0.5
drying_line_color <- fill_blue
drying_line_type <- "solid"
#####################
## plot evaporation
#####################
print(pla_dat)
print(evap_dat)

evap_draws <- decay_chains %>% spread_draws(
                                   beta[evap_class_id],
                                   drying_time[evap_class_id]) %>%
    add_titer_metadata(pla_dat, "evap_class_id") %>%
    inner_join(evap_dat %>% distinct(temperature, humidity,
                                     .keep_all = TRUE) %>%
               select(temperature,
                      humidity,
                      initial_mass,
                      equilibrium_mass),
               by = c("temperature", "humidity"))

print(evap_draws)

chosen_draws <- sample(1:max(evap_draws$.draw), n_lines_evap)

drying_time <- evap_draws %>%
    filter(temperature == chosen_temperature &
           humidity == chosen_rh) %>%
    select(drying_time, .draw)

print(drying_time)

median_drying_time <- median(
    drying_time$drying_time)

print(median_drying_time)


evap_plot_draws <- evap_draws %>%
    filter(.draw %in% chosen_draws) %>%
    filter(temperature == chosen_temperature &
           humidity == chosen_rh) %>%
    select(-time) %>%
    crossing(times) %>%
    mutate(
        pred_mass = initial_mass - beta * time,
        mass = ifelse(pred_mass >= equilibrium_mass,
                      pred_mass,
                      equilibrium_mass))

evap_plot_dat <- evap_dat %>%
    filter(temperature == chosen_temperature &
           humidity == chosen_rh)

evap_panel <- evap_plot_draws %>%
    ggplot(aes(
        x = time,
        y = mass)) +
    geom_line(aes(group = .draw),
              alpha = line_alpha_evap,
              color = "black",
              size = 1) +
    geom_point(shape = 21,
               fill = fill_blue,
               alpha = 0.7,
               size = 2.5,
               stroke = 0.75,
               data = evap_plot_dat) +
    geom_point(
        aes(y = equilibrium_mass),
        x = 24,
        shape = 22,
        fill = fill_blue,
        alpha = 1,
        size = 5,
        stroke = 0.75,
        data = evap_plot_dat) +
    geom_vline(
        xintercept = median_drying_time,
        size = drying_line_size,
        alpha = drying_line_alpha,
        color = drying_line_color,
        linetype = drying_line_type) +
    scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
    coord_cartesian(ylim = ylim_evap,
                    xlim = xlim,
                    clip = "off",
                    expand = FALSE)


cat("extracting draws for decay rates / intercepts...\n")
int_draws <- decay_chains %>%
    spread_draws(intercept[titer_id]) %>%
    add_titer_metadata(pla_dat)

cat("extracting decay rates...\n")
decay_draws <- decay_chains %>%
    spread_draws(c(transient_decay_rate,
                   decay_rate)[experiment_id]) %>%
    inner_join(drying_time,
               by = ".draw")


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
    filter(material == "Plastic")  %>%
    filter(temperature == chosen_temperature &
           humidity == chosen_rh)


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
    filter(.draw %in% chosen_draws) %>%
    filter(temperature == chosen_temperature &
           humidity == chosen_rh)


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
        predicted_titer = 10^predicted_log_titer) %>%
    filter(predicted_titer >= ylim_titer[1])

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
        data = titer_ests_draws %>%
            filter(time <= xlim[2]),
        point_size = 6,
        size = 7,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    geom_vline(
        xintercept = median_drying_time,
        size = drying_line_size,
        alpha = drying_line_alpha,
        color = drying_line_color,
        linetype = drying_line_type) +
    scale_y_log10_mathformat() +
    scale_fill_virus() +
    scale_color_virus() +
    shape_scale +
    coord_cartesian(
        xlim = xlim,
        ylim = ylim_titer,
        clip = "off",
        expand = FALSE)


## styling and assembly of panels

evap_panel <- evap_panel +
    ylab("mass (mg)") +
    xlab("") +
    theme_project() +
    labs(subtitle = paste0(
             sprintf("%s \u00B0C, %s",
                     chosen_temperature,
                     chosen_rh),
             " % relative humidity")) +
    theme(plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_blank()) +
    theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))


regression_panel <- regression_panel +
    annotate("text", x = 2.5, y = 100 * ylim_titer[2],
             label = "evaporation\nphase",
             hjust = 0.5,
             vjust = 0.5,
             size = 7,
             family = theme_family) +
    annotate("text", x = 25, y = 100 * ylim_titer[2],
             label = "quasi-equilibrium\nphase",
             hjust = 0.5,
             vjust = 0.5,
             size = 7,
             family = theme_family) +
    annotate("segment",
             x = median_drying_time,
             xend = median_drying_time,
             y = 10 * ylim_titer[2],
             yend = 1000 * ylim_titer[2],
             size = drying_line_size,
             alpha = drying_line_alpha,
             color = drying_line_color,
             linetype = drying_line_type) +
    xlab(titer_xlab) +
    ylab(titer_ylab) +
    theme_project() +
    theme(legend.position = "none") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

aligned_plots <- plot_grid(
    evap_panel,
    regression_panel,
    nrow = 2,
    align = "hv",
    axis = "tlbr",
    greedy = FALSE)

blank_plot <- ggplot() +
    geom_blank() +
    theme_nothing()

top_row <- plot_grid(
    blank_plot,
    aligned_plots,
    blank_plot,
    nrow = 1,
    rel_widths = c(0.15, 0.7, 0.15))


##################################################
## calculate posterior draws for half lives
##################################################
## get needed draws and add human readable names

cat("reading chains...\n")
draws <- chains %>%
    spread_draws(decay_rate[experiment_id],
                 transient_decay_rate[experiment_id]) %>%
    add_titer_metadata(dat, "experiment_id") %>%
    mutate(half_life_equilibrium = to_half_life(decay_rate),
           half_life_transient = to_half_life(transient_decay_rate)) %>%
    pivot_longer(c(half_life_equilibrium,
                   half_life_transient),
                 names_to = "period",
                 names_prefix = "half_life_",
                 values_to = "half_life") %>%
    mutate(period =
               factor(period,
                      levels = c(
                          "equilibrium",
                          "transient")))


cat("\n\nhalf lifes:")
print(draws %>%
      group_by(virus, temperature, humidity, period) %>%
      summarise(med = quantile(half_life, 0.5),
                q025 = quantile(half_life, 0.025),
                q975 = quantile(half_life, 0.975)))


###################################
## plot half-lives as eyeplots
###################################

temps <- unique(draws$temperature) %>% sort(decreasing = TRUE)
hums <- unique(draws$humidity) %>% sort()

hl_draws <- draws %>%
    mutate(
        humidity = ifelse(
            period == "transient",
            120,
            humidity)
    )

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
    facet_grid(cols = vars(period),
               scales = "free_x",
               space = "free_x")

## styling: no facet labels because is background plot

halflives <- halflives +
    theme_project() +
    xlab("") +
    ylab("half-life (hours)") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit("0.5", "in"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.key.size = unit(2, "lines"),
        legend.position = c(0.2, 0.85),
        legend.margin = margin(t = -10,
                               l = 10,
                               r = 10,
                               b = 10),
        legend.title = element_blank(),
        legend.box.background = element_rect(color = "black",
                                             size = 2))
####################################
## import schematic and add text
####################################
cat("Assembling schematic...\n")

hot_anchor <- 0.35
cold_anchor <- 0.475

evap_x <- 45
evap_y <- 1.3
eq_droplet_y <- 0.9

efflor_x <- 25
small_x <- 60
large_x <- 85

evap_loc <- c(evap_x, evap_y)
efflor_loc <- c(efflor_x, eq_droplet_y)
small_loc <- c(small_x, eq_droplet_y)
large_loc <- c(large_x, eq_droplet_y)


process_text_delta_x <- 12.5
process_text_y <- 1.275
origin_pad <- c(12, 0.1)

y_stretch <- 1.5

schematic_text <- ggplot(
    data = tibble(xs = c(0, 100),
                  ys = c(0, 1)),
    mapping = aes(x = xs,
                  y = ys)) +
    geom_text(
        x = evap_x - process_text_delta_x,
        y = process_text_y,
        label = "efflorescence",
        hjust = 1,
        size = 10,
        color = "black",
        family = theme_family) +
    geom_text(
        x = evap_x + process_text_delta_x + 2,
        y = process_text_y,
        label = "concentration",
        size = 10,
        hjust = 0,
        color = "black",
        family = theme_family) +
    geom_text(
        x = 52.5,
        y = 0.625,
        label = "ERH",
        size = 10,
        color = "black",
        family = theme_family) +
    schematic_curve(
        y_anchor = hot_anchor,
        mapping = aes(color = temperature),
        lineend = "round",
        ncp = 10,
        size = 5,
        data = tibble(temperature = 27)) +
    geom_text(
        aes(color = temperature),
        x = 10,
        y = hot_anchor - 0.07,
        label = "hotter",
        size = 10,
        data = tibble(temperature = 27),
        family = theme_family) +
    schematic_curve(y_anchor = cold_anchor,
                    mapping = aes(color = temperature),
                    lineend = "round",
                    ncp = 10,
                    data = tibble(temperature = 10),
                    size = 5) +
    geom_text(
          aes(color = temperature),
        x = 10,
        y = cold_anchor + 0.05,
        label = "colder",
        size = 10,
        data = tibble(temperature = 10),
        family = theme_family) +
    geom_line(
        data = tibble(
            xs = c(45, 45),
            ys = c(0, 1.1)),
        color = "black",
        linetype = "11",
        size = 5) +
    schematic_arrow(evap_loc,
                    efflor_loc,
                    origin_pad,
                    c(0, 0.1),
                    curvature = 0.25,
                    size = 1.25) +
    schematic_arrow(evap_loc,
                    large_loc,
                    origin_pad,
                    c(5, 0.2),
                    curvature = -0.25,
                    size = 1.25) +
    schematic_arrow(evap_loc,
                    small_loc,
                    origin_pad,
                    c(1, 0.1),
                    curvature = -0.1,
                    size = 1.25) +
    scale_y_continuous(breaks = c()) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(xlim = c(0, 102),
                    ylim = c(0, y_stretch),
                    expand = 0) +
    scale_color_temperature(guide = FALSE) +
    xlab("relative humidity (%)") +
    ylab("") +
    theme_project() +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(hjust = -0.5))

im_width <- 0.27
im_height <- 0.27

schematic <- ggdraw() +
    draw_plot(schematic_text) +
    draw_image(evap_img,
               x = 0.49,
               y = evap_y / y_stretch,
               width = im_width,
               height = im_height,
               hjust = 0.5,
               vjust = 0.5) +
    draw_image(efflor_img,
               x = (efflor_x + 5) / 100,
               y = (eq_droplet_y + 0.01) / y_stretch,
               width = im_width,
               height = im_height,
               hjust = 0.5,
               vjust = 0.5) +
    draw_image(small_img,
               x = small_x / 100,
               y = eq_droplet_y / y_stretch,
               width = im_width,
               height = im_height,
               hjust = 0.5,
               vjust = 0.5) +
    draw_image(large_img,
               x = large_x / 100,
               y = (eq_droplet_y + 0.045) / y_stretch,
               width = im_width,
               height = im_height,
               hjust = 0.5,
               vjust = 0.5)

####################################
## compose full figure
####################################

cat("Composing figure from subpanels...\n")


bottom_row <- plot_grid(
    halflives,
    schematic,
    ncol = 2,
    labels = c("b", "c"),
    label_size = 30)


full_plot <- plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    labels = c("a", ""),
    label_size = 30)

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')

save_plot(
    outpath = outpath,
    fig = full_plot,
    base_height = 20,
    base_asp = 0.8,
    limitsize = FALSE)

warnings()
