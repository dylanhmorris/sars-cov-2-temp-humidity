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

if( length(args) > 2) {
    data_path <- args[1]
    check_results_path <- args[2]
    outpath <- args[3]
} else {
    stop(paste0("Must at least supply data path, ",
                "check results path, ",
                "and output path"))
}

## read data / style files

cat("reading data (this may take a while)...\n")
check_chains <- readRDS(check_results_path)
dat <- read_csv(data_path,
                col_types = cols()) %>%
    mutate(

        display_material = ifelse(
            material == "Bulk medium",
            paste0("Bulk ", medium_group),
            material),
        
        name_no_hum = paste0(
            study_author_1, "\n",
            virus, "\n",
            display_material, "\n",
            temperature,
            "\u00B0C"),

        experiment_display_name =
            ifelse(!is.na(humidity),
                       paste0(name_no_hum,
                              "/",
                              humidity,
                              " %RH"),
                   name_no_hum)) %>%
    arrange(study_author_1,
            virus,
            material,
            medium_group,
            temperature) %>%
    mutate(
        experiment_display_name = factor(
            experiment_display_name,
            unique(experiment_display_name))) ## reorder


dat <- dat %>%
    group_by(experiment_id) %>%
    summarise(max_time = max(elapsed_time)) %>%
    inner_join(dat,
               by = "experiment_id")

cat("data read succesfully!\n")


#################################
## overall plot styling
#################################
set.seed(34634) # reproducible! (since we use random draws)

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")
material_colors <- get_params("material_colors")

text_size = 40
titer_ylab <- expression("fraction remaining viable")
titer_xlab <- "time (hours)"
line_alpha <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

n_lines <- 10

##################################################
## calculate posterior draws for regression lines
##################################################

tidy_draws <- check_chains %>%
    spread_draws(sampled_titers_pred[measurement_id]) %>%
    add_titer_metadata(dat, "measurement_id") %>%
    ungroup()

panel <- tidy_draws %>%
    ggplot(aes(x = elapsed_time,
               y = 10^(sampled_titers_pred),
               fill = virus)) +
    geom_hline(
        aes(yintercept = 10^log10_LOD_fraction),
        linetype = "dashed",
        data = dat %>%
            distinct(
                experiment_id,
                .keep_all = TRUE)) +
    geom_violin(
        aes(group = elapsed_time,
            width = max_time),
        position = "identity",
        alpha = 0.5) + 
    geom_point(aes(y = 10^log10_fraction_viable,
                   shape = virus),
               alpha = 0.9,
               data = dat) +
    scale_fill_virus() +
    scale_shape_virus() +
    scale_y_log10_mathformat() +
    coord_cartesian(ylim = c(10^(-6), 10)) +
    facet_wrap(vars(experiment_display_name),
               scales = "free_x",
               ncol = 5)

panel <- panel +
    theme_project(base_size = 10) +
    xlab(titer_xlab) +
    ylab(titer_ylab)

####################################
## compose full figure from panels
####################################

labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.tag=element_text(angle=-90,
                          size = text_size),
    plot.tag.position=c(1.05, 0.5))

left_margin <- theme(
    plot.margin = margin(b = 3, t = 1, l = 1, r = 0, unit = "cm"))
    
cat('making full figure...\n')

cat("n expected panels: ", dat %>%
                           select(experiment_id) %>%
                           unique() %>%
                           count() %>%
                           as.numeric(), "\n")

full_fig <- panel + labeled_panel_theme + left_margin

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 14,
          base_asp = 0.8)
warnings()
