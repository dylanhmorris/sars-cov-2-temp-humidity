#!/usr/bin/env Rscript

########################################
## filename: figure-compare-activation-energies.R
## author: Dylan Morris <dhmorris@princeton.edu>
## description: create figure showing posterior
## of difference between E_a_sol and E_a_eff
#######################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidyr',      # for pivot_longer()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # plotting
    'cowplot',    # publication-ready figures
    'extrafont',  # CM font
    'virusenv'    # project functions
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

#################################
## read command line args
#################################
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
dat_path <- args[1]
evap_dat_path <- args[2]
meas_results_path <- args[3]
mod_results_path <- args[4]
outpath <- args[5]


#################################
# set styling
#################################

stroke <- 2
interval_size_range <- c(2, 4)
alpha <- 0.9
parameter_colors <- get_params("parameter_colors")
fill_color_E_a <- parameter_colors["fill_color_E_a"]
fill_color_diff <- "#97afc2"
stackratio <- 1
param_levels <- c(
    "E_a_dry",
    "E_a_wet",
    "E_a_diff",
    "E_a_percent")

param_labels <- c("E[a]^{eff}",
                  "E[a]^{sol}",
                  "E[a]^{eff} - E[a]^{sol}",
                  "frac(E[a]^{eff} - E[a]^{sol}, E[a]^{eff})")

model_levels <- c("main model",
                  "directly measured concentration")

model_labels <- c("main~model",
                  "directly~measured~concentration")

#################################
# read in needed data
#################################

dat <- read_data_for_plotting(dat_path)

meas_results_chains <- readRDS(meas_results_path)
mod_results_chains <- readRDS(mod_results_path)


parse_chains <- function(chains){
    result <- chains %>%
        spread_draws(log_k_30_20_ratio_dry,
                     log_k_30_20_ratio_wet,
                     E_a_dry,
                     E_a_wet) %>%
        mutate(log_ratio_ratio =
                   log_k_30_20_ratio_dry - log_k_30_20_ratio_wet,
               E_a_diff = E_a_dry - E_a_wet,
               E_a_percent = E_a_diff / E_a_dry)
    
    return(result)
}

meas_results_draws <- meas_results_chains %>%
    parse_chains() %>%
    mutate(model = "directly measured concentration")

mod_results_draws <- mod_results_chains %>%
    parse_chains() %>%
    mutate(model = "main model")

results_draws <- bind_rows(meas_results_draws,
                           mod_results_draws) %>%
    mutate(model = factor(model,
                          levels = model_levels,
                          labels = model_labels))

results_melt <- results_draws %>%
    pivot_longer(c(E_a_dry, E_a_wet, E_a_diff, E_a_percent),
                 names_to = "parameter",
                 values_to = "value") %>%
    mutate(param = factor(parameter,
                          levels = param_levels,
                          labels = param_labels))

figure_E_a <- results_melt %>%
    filter(parameter %in% c("E_a_dry",
                            "E_a_wet")) %>%
    ggplot(aes(
        x = value)) +
    stat_dotsinterval(
        quantiles = 100,
        fill = fill_color_E_a,
        stroke = stroke,
        alpha = alpha,
        stackratio = stackratio,
        interval_size_range = interval_size_range,
        binwidth = 1500,
        slab_colour = "black") +
    facet_grid(vars(param),
               vars(model),
               labeller = label_parsed) +
    ylab("posterior frequency") +
    xlab("activation energy (J/mol)") +
    theme_project()

max_x <- max(c(abs(quantile(results_draws$E_a_percent, 0.003)),
               abs(quantile(results_draws$E_a_percent, 0.997))))

figure_diff <- results_melt %>%
    filter(parameter == "E_a_percent") %>%
    ggplot(aes(
        x = value)) +
    stat_dotsinterval(
        quantiles = 100,
        fill = fill_color_diff,
        stroke = stroke,
        alpha = alpha,
        binwidth = 0.033,
        stackratio = stackratio,
        interval_size_range = interval_size_range,
        slab_colour = "black") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    facet_grid(vars(param),
               vars(model),
               labeller = label_parsed) +
    ylab("posterior frequency") +
    xlab("percent difference in activation energies") +
    coord_cartesian(xlim = c(-max_x, max_x)) +
    theme_project() +
    theme(strip.text.x = element_blank())

figure <- plot_grid(
    figure_E_a,
    figure_diff,
    nrow = 2,
    rel_heights = c(2, 1.25))
    

save_plot(outpath = outpath,
          fig = figure,
          base_height = 15,
          base_asp = 1)

warnings()
