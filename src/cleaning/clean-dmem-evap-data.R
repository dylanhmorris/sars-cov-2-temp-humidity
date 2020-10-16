#!/usr/bin/env Rscript

#####################################
## name: clean_dmem_evap_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## process raw data and save cleaned
## data for evaporation of DMEM
##
####################################


suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

## functions


clean_data <- function(dat_path,
                       endpoint_path,
                       outpath) {

    cat("Reading data...\n")
    dat <- read_csv(dat_path,
                    col_names = c(
                        "time_min",
                        "T10_RH40",
                        "T10_RH65",
                        "T10_RH85",
                        "T22_RH40",
                        "T22_RH65",
                        "T22_RH85",
                        "T27_RH40",
                        "T27_RH65",
                        "T27_RH85"),
                    skip = 1)

    end_dat <- read_csv(endpoint_path) %>%
        select(-drying_time)
    
    cat("Processing data...\n")

    ## tidy data
    dat <- dat %>%
        pivot_longer(
            -time_min,
            names_to = "condition",
            values_to = "measured_mass") %>%
        mutate(time = time_min / 60,
               temperature = as.numeric(substr(condition, 2, 3)),
               humidity = as.numeric(substr(condition, 7, 8)))

    cleaned <- dat %>%
        select(time, temperature, humidity, measured_mass) %>%
        filter(!is.na(measured_mass)) %>%
        inner_join(end_dat,
                   by = c("temperature", "humidity"))
    
    cat("Writing results to ", outpath, "...\n")
    write_csv(cleaned,
              outpath)

}

## load and clean data
args <- commandArgs(trailingOnly = TRUE)
raw_data_path <- args[1]
raw_end_path <- args[2]
outpath <- args[3]

clean_data(raw_data_path,
           raw_end_path,
           outpath)

warnings()
