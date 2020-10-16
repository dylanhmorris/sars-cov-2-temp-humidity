#!/usr/bin/env Rscript

########################################
## filename: table-halflives.R
## author: Dylan Morris <dhmorris@princeton.edu>
## create table of half-life posterior medians
## and 95% CI
#######################################

script_packages <- c(
    'virusenv',   # project functions
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'xtable',     # to create LaTeX table
    'stringr'     # string regex
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}


#################################
# functions
#################################
sanitize_T <- function(temp_numeric){
    return( case_when(
        temp_numeric == 10 ~ "LowC",
        temp_numeric == 22 ~ "MidC",
        temp_numeric == 27 ~ "HighC"))
}

sanitize_RH <- function(RH_numeric){
    return( case_when(
        RH_numeric == 40 ~ "LowRH",
        RH_numeric == 65 ~ "MidRH",
        RH_numeric == 85 ~ "HighRH"))
}


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
dat_path <- args[1]
hl_path <- args[n_args - 1] 
outpath <- args[n_args]

dat <- read_data_for_plotting(dat_path)

save_table <- grepl("table", outpath)
save_macros <- grepl("macros", outpath)


hl_chains <- readRDS(hl_path)

hl_draws <- hl_chains %>%
    spread_draws(log_half_life[experiment_id],
                 transient_decay_rate[experiment_id]) %>%
    mutate(transient_log_half_life = transient_decay_rate %>%
               to_half_life() %>%
               log()) %>%
    add_titer_metadata(dat, "experiment_id")

macro_list <- list()
macro_template <- "\\newcommand{\\%s}{%s}\n"

model_name <- str_match(outpath,
                        "(macros-|table-)([aA-zZ|-]+)-halflives.*")[3]

if(save_table) {
    cat("Half-life table for model",
        model_name,
        "...\n")
} else if(save_table){
    cat("Half-life macros for model",
        model_name,
        "...\n")
}

if(model_name == "sars-mers"){
    viruses <- "SARS-CoV-1 and MERS-CoV"
} else {
    viruses <- "SARS-CoV-2"
}
                          

    
##################################################
## calculate half-life table
##################################################

table_caption <- paste0(
    sprintf(
        paste0(
            "\\textbf{Estimated half-lives in hours of %s ",
            "on polypropylene as a function of temperature (T) ",
            "and relative humidity (RH).} "),
        viruses),
    
    "Estimated half-lives are reported as posterior ",
    "median and the middle 95\\% credible interval.")

table_label <- sprintf("tab:%s-halflives",
                       model_name)
                       

hl_draws <- hl_draws %>%
    mutate(hl = exp(log_half_life),
           thl = exp(transient_log_half_life))

hl_table <- hl_draws %>%
    group_by(temperature, humidity, virus) %>%
    summarise(median = median(hl),
              q025 = quantile(hl, 0.025),
              q975 = quantile(hl, 0.975),
              ) %>%
    arrange(temperature, humidity) %>%
    mutate(humidity = factor(humidity))


if(model_name != "sars-mers")
    hl_table <- hl_table %>%
        select(-virus)

n_columns <- dim(hl_table)[2]
n_word_columns <- n_columns - 3
table_precisions <- c(rep(0, n_word_columns + 1),
                      rep(2, 3))


for(row in list(hl_table)){
    
    med_macro_name <- paste0("EqHalfLifeMedian",
                             sanitize_T(row$temperature),
                             sanitize_RH(row$humidity))
    q025_macro_name <- paste0("EqHalfLifeLowQuant",
                              sanitize_T(row$temperature),
                              sanitize_RH(row$humidity))
    q975_macro_name <- paste0("EqHalfLifeHighQuant",
                              sanitize_T(row$temperature),
                              sanitize_RH(row$humidity))
    macro_list[med_macro_name] <- row$median
    macro_list[q025_macro_name] <- row$q025
    macro_list[q975_macro_name] <- row$q975
}

transient_table <- hl_draws %>%
    filter(humidity == 40) %>%
    group_by(temperature, virus) %>%
    summarise(median = median(thl),
              q025 = quantile(thl, 0.025),
              q975 = quantile(thl, 0.975),
              ) %>%
    ungroup() %>%
    mutate(humidity = factor("")) %>%
    arrange(temperature)


if(model_name != "sars-mers")
    transient_table <- transient_table %>%
        select(-virus)


for(row in list(transient_table)){
    
    med_macro_name <- paste0("EvaporationHalfLifeMedian",
                             sanitize_T(row$temperature))
    q025_macro_name <- paste0("EvaporationHalfLifeLowQuant",
                              sanitize_T(row$temperature))
    q975_macro_name <- paste0("EvaporationHalfLifeHighQuant",
                              sanitize_T(row$temperature))
    macro_list[med_macro_name] <- row$median
    macro_list[q025_macro_name] <- row$q025
    macro_list[q975_macro_name] <- row$q975
}


full_table <- hl_table %>%
    bind_rows(transient_table)


print(full_table)

num_eq <- dim(hl_table)[1]

rowformatter <- function(rownames){
    result <- ifelse(rownames == 1,
                     "quasi-equilibrium phase",
              ifelse(rownames == num_eq + 1,
                     "evaporation phase",
                     ""))

    return( result )
}



latex_table <- full_table %>%
    rename("T (\\si{\\celsius})" = temperature,
           "RH (\\si{\\%})" = humidity,
           "median half-life (h)" = median,
           "2.5 \\si{\\%}" = q025,
           "97.5 \\si{\\%}" = q975) %>%
    xtable(
        digits = table_precisions,
        caption = table_caption,
        label = table_label)

####################################
## save table / macros
####################################

if(save_table){
    cat('saving table to', outpath, '...\n')

    if(model_name == "sars-mers"){
        rowbreaks <- list()
        rowbreaks$pos <- list(c(2))
        rowbreaks$command <- "\\\\ \n"
    } else {
        rowbreaks <- list()
        rowbreaks$pos <- list(c(3, 6, 9))
        rowbreaks$command <- "\\\\ \n"
    }
    
    print(latex_table,
          caption.placement = "top",
          size = "\\small",
          sanitize.text.function = function(x){x},
          sanitize.rownames.function = rowformatter,
          add.to.row = rowbreaks,
          math.style.exponents = TRUE,
          file = outpath)
} else if (save_macros) {

    cat('saving macros to', outpath, '...\n')

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
