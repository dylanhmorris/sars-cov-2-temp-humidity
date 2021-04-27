#!/usr/bin/env Rscript

########################################
## filename: figure-predictive-check.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot predictive checks for T/RH stability
#######################################

suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(extrafont))


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)

n_args <- length(args)

outpath <- args[n_args]

##
## determine what is being computed
##
prior <- grepl("prior", outpath)
posterior <- grepl("posterior", outpath)
infer <- grepl("inferred", outpath)
evaporation <- grepl("evap-phase", outpath)

if( prior ) {
    check_results_path <- args[n_args - 1]
    titers_path <- args[n_args - 2]
    evap_path <- args[n_args - 3]
    data_path <- args[1]
} else if( posterior ){
    check_results_path <- args[n_args - 1]
    titers_path <- args[n_args - 2]
    evap_path <- check_results_path
    data_path <- args[1]
}

## read data / style files

cat("reading data (this may take a while)...\n")
check_chains <- readRDS(check_results_path)
ests_chains <- readRDS(titers_path)
evap_chains <- readRDS(evap_path)

dat <- read_csv(data_path,
                col_types = cols())

cat("data read succesfully!\n")


#################################
## overall plot styling
#################################
set.seed(34634) # reproducible! (since we use random draws)

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")
material_colors <- get_params("material_colors")

text_size = 40
detection_linesize = 0.75
ylim <- c(1e-3, 1e6)
LOD_log10_per_ml <- 0.5
LOD <- 10^LOD_log10_per_ml
line_alpha <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

conversion_factor <- 1
n_lines <- 10

if(evaporation | infer){
    xlab_text <- "time since deposition (hours)"
} else {
    xlab_text <- "time since quasi-equilibrium reached (hours)"
}

print(xlab_text)

if(grepl("sars-mers", outpath)){
    
    faceting <- function() {
        return( list(facet_wrap(vars(virus),
                                ncol = 2)))
    }
    
    plot_labels <- function(){
        return(
            list(
                labs(subtitle = "virus"),
                xlab(xlab_text),
                ylab(expression("virus titer (TCID"[50] * "/mL media)")))
        )
    }

} else {
    faceting <- function() {
        return( list(
                   facet_grid(vars(humidity),
                              vars(temperature))))
    }
    
    plot_labels <- function(){
        return(
            list(
                labs(tag = "relative humidity (%)",
                     subtitle = "temperature (\u00B0C)"),
                xlab(xlab_text),
                ylab(expression("virus titer (TCID"[50] * "/mL media)")))
        )
    }
}


## calculate posterior median drying times

median_drying_times <- evap_chains %>%
    spread_draws(drying_time[evap_class_id]) %>%
    group_by(evap_class_id) %>%
    summarise(med_drying_time = median(drying_time))

dat <- dat %>%
    inner_join(median_drying_times,
               by = "evap_class_id")

## get needed draws and add human readable names

if(evaporation){
    xlim <- c(0, max(dat$med_drying_time))
} else {
    xlim <- c(0, max(dat$time - dat$med_drying_time))
}

print(xlim)

autoscale <- xlim[2] / 2

##################################################
## calculate posterior draws for regression lines
##################################################

cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")


tidy_draws <- check_chains %>%
    spread_draws(sampled_titers_pred[titer_id]) %>%
    add_titer_metadata(dat) %>%
    ungroup()

cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")

if(evaporation){
    tidy_draws <- tidy_draws %>%
        filter(time < med_drying_time)
} else {
    tidy_draws <- tidy_draws %>%
        filter(time >= med_drying_time) %>%
        mutate(time = time - med_drying_time)
}

cat("extracting positive wells...\n")
pos_wells <- dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

cat('extracting titer estimates...\n')

titer_ests_draws <- ests_chains %>%
    spread_draws(sampled_titer[titer_id])

## get human readable names and detectability
titer_ests_draws <- titer_ests_draws %>%
    add_titer_metadata(dat) %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 0) %>%
    filter(material == "Plastic")

## filter time
if(evaporation){
    print("evaporation")
    titer_ests_draws <- titer_ests_draws %>%
        filter(time < med_drying_time) %>%
        arrange(desc(time))
    tidy_draws <- tidy_draws %>%
        filter(time < med_drying_time)
} else {
    titer_ests_draws <- titer_ests_draws %>%
        filter(time >= med_drying_time) %>%
        mutate(time = time - med_drying_time) %>%
        arrange(desc(time))
}

titer_ests_draws <- titer_ests_draws %>%
    mutate(
        log10_titer_per_ml = ifelse(
            detectable,
            sampled_titer + 1,
            LOD_log10_per_ml))



## adding one to convert to per mL from per 0.1 mL

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

panel <- tidy_draws %>%
    ggplot(aes(x = time,
               y = 10^(1 + sampled_titers_pred),
               fill = virus)) +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    geom_violin(
        aes(group = time),
        width = autoscale,
        position = "identity",
        alpha = 0.5) + 
    stat_pointinterval(
        mapping = aes(x = time,
                      y = 10^log10_titer_per_ml,
                      shape = detectable,
                      fill = virus,
                      group = titer_id),
        .width = 0.95,
        fatten_point = 2,
        interval_size = 14,
        stroke = 1,
        data = titer_ests_draws) +
    scale_fill_virus() +
    scale_fill_virus(aesthetics = "point_fill") +
    shape_scale + 
    scale_y_log10_mathformat() +
    coord_cartesian(xlim = xlim,
                    ylim = ylim) +
    faceting() +
    plot_labels() +
    theme_project(base_size = text_size) +
    theme(legend.position = "none",
          panel.spacing.y = unit(5, "lines"))

####################################
## compose full figure from panels
####################################

labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.tag = element_text(angle = -90,
                            size = text_size),
    plot.tag.position=c(1.05, 0.5))

left_margin <- theme(
    plot.margin = margin(b = 3, t = 1, l = 1, r = 4, unit = "cm"))
    
cat('making full figure...\n')

full_fig <- panel + labeled_panel_theme + left_margin

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 18,
          base_asp = 1.2)
warnings()
