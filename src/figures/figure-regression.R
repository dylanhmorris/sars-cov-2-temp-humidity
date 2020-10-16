#!/usr/bin/env Rscript

########################################
## filename: figure-regression.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot main regression figure for environmental
## conditions analysis
#######################################

suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(extrafont))


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
evap_path <- args[2]
results_path <- args[3]
titers_path <- args[4]
outpath <- args[5]

## read data / style files

cat("reading data (this may take a while)...\n")
decay_chains <- readRDS(results_path)
ests_chains <- readRDS(titers_path)

dat <- read_data_for_plotting(data_path)
evap_dat <- read_csv(evap_path,
                     col_types = cols())
cat("data read successfully!\n")

## model to use
evaporation_with_concentration <- (
    any(grepl("fraction_solute", names(decay_chains))))

evaporation_biphasic <- (
    (!evaporation_with_concentration) &
    (any(grepl("transient_decay_rate", names(decay_chains)))))

#################################
## overall plot styling
#################################

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")

set.seed(23) # reproducible! (since we use random draws)
text_size = 40
detection_linesize = 0.75
titer_ylab <- expression("virus titer (TCID"[50] * "/mL media)")
ylim <- c(1, 5e4)
LOD_log10 <- 0.5
LOD <- 10^LOD_log10

conversion_factor <- 1
n_lines <- 10
line_alpha <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

evaporation <- grepl("evap-phase", outpath)

##################################################
## calculate posterior draws for regression lines
##################################################

## calculate posterior median drying times

median_drying_times <- decay_chains %>%
    spread_draws(drying_time[experiment_id]) %>%
    group_by(experiment_id) %>%
    summarise(med_drying_time = median(drying_time))

## get needed draws and add human readable names

if(evaporation){
    xlim <- c(0, max(median_drying_times$med_drying_time))
} else {
    xlim <- c(0, 96)
}

fineness = 250
plot_times <- tibble(
    time = seq(xlim[1], xlim[2],
               length.out = fineness))


cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")
int_draws <- decay_chains %>%
    spread_draws(intercept[titer_id]) %>%
    add_titer_metadata(dat)


if( evaporation_with_concentration ){
    cat("Using explicit concentration model...\n")
    decay_draws <- decay_chains %>%
        spread_draws(c(transient_decay_rate,
                       decay_rate)[experiment_id],
                     initial_mass_fraction_solute)

    evap_draws <- decay_chains %>%
        spread_draws(c(beta,
                       drying_time)[evap_class_id])

    predict_log_titers <- function(data){
        return( data %>%
                mutate(predicted_log_titer =
                           1 + predict_titers_explicit_evaporation(
                                   time,
                                   intercept,
                                   transient_decay_rate,
                                   decay_rate,
                                   beta / (initial_mass * (1 - initial_mass_fraction_solute)),
                                   equilibrium_concentration_factor)))
    }

} else if( evaporation_biphasic ) {
    cat("Using biphasic evaporation model...\n")
    decay_draws <- decay_chains %>%
        spread_draws(c(transient_decay_rate,
                       decay_rate)[experiment_id])
    evap_draws <- decay_chains %>%
        spread_draws(c(beta,
                       drying_time)[evap_class_id]) %>%
        ungroup()

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
}

cat("extracting positive wells...\n")
pos_wells <- dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

tidy_draws <- decay_draws %>%
    inner_join(int_draws %>%
               select(intercept,
                      experiment_id,
                      .draw,
                      evap_class_id,
                      temperature,
                      humidity,
                      titer_id,
                      virus),
               by = c("experiment_id",
                      ".draw")) %>%
    inner_join(evap_draws %>%
               select(evap_class_id,
                      drying_time,
                      beta,
                      .draw),
               by = c("evap_class_id",
                      ".draw")) %>%
    inner_join(evap_dat %>% distinct(temperature, humidity,
                                     .keep_all = TRUE) %>%
               select(-time),
               by = c("temperature", "humidity"))

print(tidy_draws)

if(evaporation_with_concentration){
    tidy_draws <- tidy_draws %>%
        mutate(
            equilibrium_concentration_factor = exp(mass_change_to_log_concentration_factor(
                equilibrium_mass,
                initial_mass,
                initial_mass_fraction_solute)))

}

cat('extracting titer estimates...\n')

titer_ests_draws <- ests_chains %>%
    spread_draws(sampled_titer[titer_id])

## get human readable names and detectability
titer_ests_draws <- titer_ests_draws %>%
    add_titer_metadata(dat) %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 1) %>%
    filter(material == "Plastic") %>%
    inner_join(median_drying_times,
               by = "experiment_id")

## filter time
if(evaporation){
    print("evaporation")
    titer_ests_draws <- titer_ests_draws %>%
        filter(time < med_drying_time) %>%
        mutate(time_use = time) %>%
        arrange(desc(time_use))
    figure_xlab <- "time since deposition (hours)"
} else {
    titer_ests_draws <- titer_ests_draws %>%
        filter(time >= med_drying_time) %>%
        mutate(time_use = time - med_drying_time) %>%
        arrange(desc(time_use))
    figure_xlab <- "time since quasi-equilibrium reached (hours)"
}

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
chosen_draws <- sample(1:max(tidy_draws$.draw), n_lines)

func_samples <- tidy_draws %>%
    filter(.draw %in% chosen_draws)

## annotate lines so that each
## has a unique id for ggplot overplotting
## (else two lines from the same draw but
## different replicates can get confused
## with each other)

## cross product decay_rates with x (time) values
## and calculate y (titer) values
cat('setting up x values...\n')

to_plot <- func_samples %>%
    crossing(plot_times)

if(evaporation){
    to_plot <- to_plot %>%
        filter(time < drying_time) %>%
        mutate(time_use = time)
} else {
    to_plot <- to_plot %>%
        filter(time >= drying_time) %>%
        mutate(time_use = time - drying_time)
}

## adding one to convert to per mL from per 0.1 mL
to_plot <- to_plot %>%
    predict_log_titers() %>%
    mutate(
        predicted_titer = 10^predicted_log_titer) %>%
    filter(predicted_titer > ylim[1])

shape_scale <- scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))


hl_dat <- dat %>%
    distinct(experiment_id,
             .keep_all = TRUE) %>%
    select(experiment_id,
           material,
           temperature,
           virus,
           humidity) %>% 
    inner_join(decay_draws,
               by = "experiment_id") %>%
    mutate(
        half_life = log10(2) / decay_rate,
        tenfold = log10(10) / decay_rate)


plot_dat <- titer_ests_draws

panel <- to_plot %>%
    ggplot() +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    geom_line(aes(
        x = time_use,
        y = predicted_titer,
        color = virus,
        group = interaction(.draw, titer_id)),
        alpha = line_alpha) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(x = time_use,
                      y = 10^log10_tcid50,
                      shape = detectable,
                      fill = virus,
                      group = titer_id),
        data = plot_dat,
        point_size = 6,
        size = 7,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9)

panel <- panel +
    scale_fill_virus() +
    scale_fill_virus(aesthetics = "point_fill") +
    scale_color_virus() +
    shape_scale + 
    scale_y_log10_mathformat(expand = c(0, 0)) +
    coord_cartesian(ylim = ylim,
                    xlim = xlim,
                    clip = "off") +
    facet_grid(rows = vars(humidity),
               cols = vars(temperature))


## styling: no facet labels because is background plot

panel <- panel +
    theme_project(base_size = text_size) +
    theme(panel.border = element_rect(size = 2,
                                      color = "black",
                                      fill = NA)) +
    xlab(figure_xlab) +
    ylab(titer_ylab) +
    theme(legend.position = "none",
          panel.spacing.y = unit(5, "lines")) +
    labs(tag = "relative humidity (%)",
         subtitle = "temperature (\u00B0C)")


if(!evaporation){
    panel <- panel +
        geom_density_ridges(
            data = hl_dat,
            mapping = aes(
                y = 10^(log10(ylim[2])),
                x = half_life,
                fill = virus,
                height = ..ndensity..),
            scale = 1.2) +
        stat_pointinterval(
            data = hl_dat,
            mapping = aes(
                y = 10^(log10(ylim[2])),
                x = half_life),
            .width = 0.95,
            shape = 19,
            color = interval_grey,
            size = 35)
}


####################################
## compose full figure from panels
####################################

labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.tag = element_text(angle=-90,
                          size = text_size),
    plot.tag.position = c(1.05, 0.5))

left_margin <- theme(
    plot.margin = margin(b = 3, t = 1, l = 1, r = 3, unit = "cm"))
    
cat('making full figure...\n')

full_fig <- panel + labeled_panel_theme + left_margin

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 15,
          base_asp = 1.2)
warnings()
