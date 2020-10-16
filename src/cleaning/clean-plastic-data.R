#!/usr/bin/env Rscript

#####################################
## name: clean_plastic_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## process raw data and save cleaned
## data for use in model fitting
## of the temp/humidity models
## for plastic
##
####################################


suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

## functions


clean_well_data <- function(dat,
                            viruses) {
    
    dat <- dat %>%
        mutate(temperature = temperature %>%
                   recode("21-23" = "22",
                          "25-27" = "27",
                          "10" = "10") %>%
                   as.character(),
               target_temperature = temperature %>%
                   as.character() %>%
                   as.numeric())
    
    dat <- dat %>% filter(material == "Plastic")
        
    ## for now, do not estimate dose deposited
    ## loss rate -- maybe add if everything
    ## else works

    dat <- dat %>% filter(!grepl("cytotoxic", comment) &
                          virus %in% viruses)

    dat$treatment <- dat$treatment %>%
        replace_na("Untreated")

    dat <- dat %>%
        mutate(caron_chamber = ((humidity == 65) |
                                (temperature == 10)))

    dat <- dat %>%
        filter(treatment == "Untreated")

    dat <- dat %>%
        group_by(virus,
                 temperature,
                 humidity,
                 material) %>%
        mutate(experiment_id = cur_group_id())
    dat <- dat %>%
        group_by(temperature,
                 humidity) %>%
        mutate(evap_class_id = cur_group_id())
    
    dat <- dat %>%
        group_by(temperature, virus, material) %>%
        arrange(temperature, virus, material, .by_group = TRUE) %>%
        mutate(transient_class_id = cur_group_id())

    dat <- dat %>%
        arrange(material) %>%
            group_by(material) %>%
        mutate(material_id = cur_group_id())
    
    dat <- dat %>%
        group_by(virus,
                 material,
                 temperature,
                 humidity,
                 replicate,
                 time) %>%
        mutate(titer_id = cur_group_id())

    dat <- dat %>%
        group_by(virus,
                 material,
                 temperature,
                 humidity,
                 time) %>%
        mutate(replicate_set_id = cur_group_id()) %>%
        ungroup()

    dat <- dat %>%
        mutate(target_humidity = humidity / 100)
    
    dat <- dat %>% arrange(titer_id)
    
    print(dat)
    return (dat)
}

## load and clean data
args <- commandArgs(trailingOnly = TRUE)
raw_data_path <- args[1]
outpath <- args[2]

delim <- ";"

cat("specifying column types...\n")

col_types <- cols_only(
    virus = col_factor(),
    material = col_factor(),
    treatment = col_character(),
    time = col_double(),
    dilution = col_integer(),
    replicate = col_integer(),
    virus_detect = col_integer(),
    temperature = col_factor(),
    humidity = col_double(),
    comment = col_character())


cat("Reading raw data...\n")
dat <- read_delim(raw_data_path,
                  delim = delim,
                  col_types = col_types)

viruses <- c("SARS-CoV-2")

if(grepl("sars-mers", outpath))
    viruses <- c("SARS-CoV-1", "MERS-CoV")

cat("cleaning well data...\n")
cleaned <- clean_well_data(dat,
                           viruses)

cat("Writing results to ", outpath, "...\n")
write_csv(cleaned,
          outpath)

warnings()
