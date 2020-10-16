#!/usr/bin/env Rscript

#####################################
## name: clean_literature_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## clean virus data from the literature
## to prepare for model fitting
####################################

suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

dat_path <- args[1]
out_path <- args[2] 

dat <- read_delim(dat_path, delim = ";")

## collect only rows of interest
dat <- dat %>%
    group_by(study_doi, virus, material, medium, humidity) %>%
    mutate(two_temps = n_distinct(temperature) > 1) %>%
    ungroup() %>%
    group_by(study_doi, virus, material, medium, temperature) %>%
    mutate(two_hums = n_distinct(humidity) > 1,
           has_hum = n_distinct(humidity) > 0) %>%
    ungroup() %>%
    filter(!grepl("half-life", estimate)) %>%
    filter(grepl("SARS", virus) |
           grepl("MERS", virus) |
           grepl("HCoV", virus)) %>%
    filter(material != "Aerosols") %>%
    filter(study_doi != "10.1155/2011/734690") %>%
    filter(medium_include) %>%
    filter(!is.na(value)) %>%
    filter(two_temps | has_hum) %>%
    rowwise() %>%
    mutate(unit_LOD = grepl(unit, LOD) | LOD == "Unspecified") %>%
    filter(unit_LOD) ## this makes sure that LODs have same units as values


## transform value and LOD to numeric columns
dat <- dat %>%
    mutate(
        LOD = as.numeric(str_extract(LOD, "[0-9]+\\.*[0-9]*")),
        value = as.numeric(ifelse(
            value == "<LOD",
            LOD,
            value)))

## create temp study ids
dat <- dat %>%
    group_by(study_doi) %>%
    mutate(study_id = cur_group_id()) %>%
    ungroup() %>%
    group_by(study_id, virus, temperature,
             humidity, material, medium) %>%
    mutate(experiment_id = cur_group_id()) %>%
    ungroup()

## determine mean at t0
dat <- dat %>%
    group_by(experiment_id) %>%
    filter(n() > 1) %>%
    mutate(elapsed_time = time - min(time)) %>%
    ungroup()

dat <- dat %>%
    filter(elapsed_time == 0) %>%
    group_by(experiment_id) %>%
    summarise(mean_zero = mean(value)) %>%
    ungroup() %>%
    right_join(dat, by = "experiment_id") %>%
    ungroup() %>%
    filter(elapsed_time > 0)


## convert everything to log10 faction viable
dat <- dat %>%
    mutate(

        log10_LOD_fraction = case_when(
            is.na(LOD) ~ -99,
            unit %in% c(
                "log10TCID50/mL",
                "log10TCID50/0.05mL",
                "log10PFU/mL") ~ LOD - mean_zero,
            unit == "-log10Nt/N0" ~ LOD,
            unit == "%Viable" ~ log10(LOD / 100)),
            
        log10_fraction_viable = case_when(
            unit %in% c(
                "log10TCID50/mL",
                "log10TCID50/0.05mL",
                "log10PFU/mL") ~ value - mean_zero,
            unit == "-log10Nt/N0" ~ -1.0 * value,
            unit == "%Viable" ~ log10(value / 100))) %>%
    mutate(log10_fraction_viable = log10(10^log10_fraction_viable))


## redo study ids
dat <- dat %>%
    group_by(study_doi) %>%
    mutate(study_id = cur_group_id()) %>%
    group_by(study_id, virus, material, medium) %>%
    mutate(experiment_class_id = cur_group_id()) %>%
    ungroup() %>%
    group_by(experiment_class_id,
             temperature,
             humidity) %>%
    mutate(experiment_id = cur_group_id()) %>%
    ungroup()


## assign measurement IDs
dat <- dat %>%
    group_by(experiment_id,
             elapsed_time) %>%
    mutate(measurement_id = cur_group_id()) %>%
    ungroup() %>%
    distinct(measurement_id, .keep_all = TRUE)

## choose needed columns

dat <- dat %>%
    select(
        study_id,
        study_author_1,
        study_doi,
        experiment_id,
        two_hums,
        two_temps,
        measurement_id,
        experiment_class_id,
        virus,
        material,
        medium,
        medium_group,
        sealed,
        temperature,
        humidity,
        elapsed_time,
        mean_zero,
        log10_fraction_viable,
        log10_LOD_fraction)

cat("Saving cleaned data to ", out_path, "...\n") 
write_csv(dat, out_path)

warnings()
