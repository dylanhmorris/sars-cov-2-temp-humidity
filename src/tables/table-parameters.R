#!/usr/bin/env Rscript

########################################
## filename: table-parameters.R
## author: Dylan Morris <dhmorris@princeton.edu>
## create table of mechanistic
## description: parameter posterior medians
## and 95% CI
#######################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidyr',      # for pivot_longer()
    'tidybayes',  # for for spread_draws(), etc.
    'xtable',     # LaTeX tables
    'stringr',    # string regex
    'virusenv'    # project functions
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

## read command line args
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
dat_path <- args[1]
evap_dat_path <- args[2]
results_path <- args[3]
outpath <- args[4]
model_name <- args[5]


#################################
# table setup
#################################

model_name_macro_string <- sprintf(
    "\\%sconcentration{}",
    model_name)

if(model_name == "modeled") {
    model_description_string <- paste0(
        "using a fitted curve relating RH to concentration ",
        "factor, as in the main text. ")
} else {
    model_description_string <- paste0(
        "using concentration ",
        "factors directly measured in evaporation ",
        "experiments. ")
}

table_caption <- paste0(
    "Parameter estimates for the mechanistic model ",
    "of SARS-CoV-2 inactivation as a function of temperature ",
    "and humidity, ",
    model_description_string,
    "Estimates are reported as posterior ",
    "median and the middle 95\\% credible interval.")

table_label <- sprintf("tab:mech-parameters-%s",
                       model_name)

table_precisions <- c(
    0,
    0,
    0,
    2,
    2,
    2)

#' as_num
#'
#' wrap vector entries within \num{}
#' for latex siunitx package
#'
#' @param vector vector whose entries
#' are to be wrapped
#'
#' @return character vector of \num{x_i}
#' for x_i in the original vector
as_num <- function(vector){
    return( paste0("\\num{", vector, "}") )
}

#################################
## macro setup
#################################
macro_list <- list()
macro_template <- "\\newcommand{\\%s}{%s}\n"

#################################
# read in needed data
#################################

dat <- read_data_for_plotting(dat_path)

results_chains <- readRDS(results_path)

results_spread <- results_chains %>%
    spread_draws(E_a_dry,
                 log_A_dry,
                 log_A_wet) %>%
    mutate(
        A_dry = exp(log_A_dry),
        A_wet = exp(log_A_wet),
        E_a = E_a_dry)

results_draws <- results_spread %>%
    pivot_longer(c(E_a, A_dry, A_wet),
                 names_to = "parameter",
                 values_to = "value")

table <- results_draws %>%
    group_by(parameter) %>%
    summarise(median = as_num(median(value)),
              q025 = as_num(quantile(value, 0.025)),
              q975 = as_num(quantile(value, 0.975)),
              ) %>%
    ungroup() %>%
    arrange(parameter)


latex_table <- table %>%
    rename("median" = median,
           "2.5 \\si{\\%}" = q025,
           "97.5 \\si{\\%}" = q975) %>%
    mutate(
        unit = case_when(
            parameter == "E_a" ~ "\\si{\\J\\per\\mol}",
            parameter == "A_dry" ~ "\\si{\\per\\hour}",
            parameter == "A_wet" ~ "\\si{\\per\\hour}"),
        parameter = case_when(
               parameter == "E_a" ~ "\\Ea{}",
               parameter == "A_dry" ~ "\\Adry{}",
               parameter == "A_wet" ~ "\\Awet{}")) %>%
    xtable(
        digits = table_precisions,
        caption = table_caption,
        label = table_label)

save_table <- grepl("table", outpath)
save_macros <- grepl("macros", outpath)

if (save_table){
    print(latex_table,
          caption.placement = "top",
          size = "\\small",
          math.style.exponents = TRUE,
          sanitize.text.function = function(x){x},
          include.rownames = FALSE,
          file = outpath)
} else if (save_macros) {
    cat('saving macros to', outpath, '...\n')

    macro_list[paste0("MedianEa", model_name)] <-
        median(results_spread$E_a)
    macro_list[paste0("QLowEa", model_name)] <-
        quantile(results_spread$E_a, 0.025)
    macro_list[paste0("QHighEa", model_name)] <-
        quantile(results_spread$E_a, 0.975)
    
    if (file.exists(outpath))
        file.remove(outpath)

    for(param_name in names(macro_list)){
        cat(sprintf(macro_template,
                    param_name,
                    macro_list[param_name]),
            file = outpath,
            append = TRUE)
    }
}

warnings()
