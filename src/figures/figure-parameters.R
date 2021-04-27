#!/usr/bin/env Rscript

########################################
## filename: figure-parameters.R
## author: Dylan Morris <dhmorris@princeton.edu>
## description: create figure showing posterior distributions
## for key mechanistic model parameters
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
interval_size_range <- c(1, 2)
alpha <- 0.9
parameter_colors <- get_params("parameter_colors")
fill_color_A <- parameter_colors["fill_color_A"] 
fill_color_E_a <- parameter_colors["fill_color_E_a"]
stackratio <- 1

model_levels <- c("modeled concentration",
                  "measured concentration")

model_labels <- c("Main~model",
                  "`Directly-measured`~concentration")

param_levels <- c("E_a",
                  "A_dry",
                  "A_wet")
param_labels <- c("E[a]~(J/mol)",
                  "A[eff]~(1/h)",
                  "A[sol]~(1/h)")

#################################
# read in needed data
#################################

dat <- read_data_for_plotting(dat_path)

meas_results_chains <- readRDS(meas_results_path)
mod_results_chains <- readRDS(mod_results_path)


parse_chains <- function(chains){
    result <- chains %>%
        spread_draws(E_a_dry,
                     log_A_dry,
                     log_A_wet) %>%
        mutate(
            A_dry = exp(log_A_dry),
            A_wet = exp(log_A_wet),
            E_a = E_a_dry) %>%
        pivot_longer(c(E_a, A_dry, A_wet),
                     names_to = "parameter",
                     values_to = "value") %>%
        mutate(
            param = parameter,
            parameter = factor(parameter,
                               levels = param_levels,
                               labels = param_labels))

    return(result)
}

meas_results_draws <- meas_results_chains %>%
    parse_chains() %>%
    mutate(model = "measured concentration")

mod_results_draws <- mod_results_chains %>%
    parse_chains() %>%
    mutate(model = "modeled concentration")

results_draws <- bind_rows(meas_results_draws,
                           mod_results_draws) %>%
    mutate(model = factor(model,
                          levels = model_levels,
                          labels = model_labels))

figure_A <- results_draws %>%
    filter(param != "E_a") %>%
    ggplot(aes(
        x = value)) +
    stat_dotsinterval(
        quantiles = 100,
        fill = fill_color_A,
        stroke = stroke,
        alpha = alpha,
        stackratio = stackratio,
        interval_size_range = interval_size_range,
        slab_colour = "black") +
    scale_x_log10_mathformat() +
    facet_grid(vars(parameter),
               vars(model),
               labeller = label_parsed) +
    xlab("") +
    ylab("posterior frequency") +
    theme_project() +
    theme(axis.title.y = element_text(
              hjust = -0.25))



figure_E_a <- results_draws %>%
    filter(param == "E_a") %>%
    ggplot(aes(
        x = value)) +
    stat_dotsinterval(
        quantiles = 100,
        fill = fill_color_E_a,
        stroke = stroke,
        alpha = alpha,
        stackratio = stackratio,
        interval_size_range = interval_size_range,
        slab_colour = "black") +
    facet_grid(vars(parameter),
               vars(model),
               labeller = label_parsed) +
    xlab("parameter value") +
    ylab("") +
    theme_project() +
    theme(strip.text.x = element_blank())



figure <- plot_grid(
    figure_A,
    figure_E_a,
    nrow = 2,
    rel_heights = c(2, 1))

save_plot(outpath = outpath,
          fig = figure,
          base_height = 15,
          base_asp = 1.25)

warnings()
