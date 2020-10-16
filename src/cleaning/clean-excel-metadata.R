#!/usr/bin/env Rscript

###################################
## filename: clean_excel_metadata.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: extract and clean metadata
## for aerosol experiments from standard
## datasheets in XLSX format
####################################

script_packages <- c(
    'readr',    # csv output
    'readxl',   # excel input
    'dplyr',    # mutate, pipe %>%, join
    'tibble',   # nice data tables
    'stringr',  # easy string regex
    'lubridate' # datetime handling
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


extract_target_T <- function(excel_file,
                             T_cell = "H53:H54",
                             T_RH_cell = NULL,
                             valid_Ts = c("10",
                                          "22",
                                          "21-23",
                                          "27")){

    if(!is.null(T_RH_cell)){
        T <- read_excel(excel_file,
                        range = T_RH_cell)
        T <- T %>%
            str_extract("[0-9]+C") %>%
            str_extract("[0-9]+")
    } else {
        T <- read_excel(excel_file,
                        range = T_cell)
        T <- T %>%
            str_extract("[0-9]+")
    }
    
    ## validate T
    T_valid <- T %in% valid_Ts

    if(!T_valid){
        cat(sprintf(
            "Invalid T\n T = %s\n", T))
        result <- NA
    } else {
        result <- T
    }
    return( result )
}
        
extract_target_RH <- function(excel_file,
                              RH_cell = "D53:D54",
                              T_RH_cell = NULL,
                              valid_RHs = c(0.4, 0.65, 0.85)){

    if(!is.null(T_RH_cell)){
        RH <- read_excel(excel_file,
                         range = T_RH_cell)
        RH <- RH %>%
            str_extract("[0-9]+%") %>%
            str_extract("[0-9]+") %>%
            as.numeric() / 100            
    } else {
        RH <- as.numeric(
            read_excel(excel_file,
                       range = RH_cell))
    }
    ## validate RH
    RH_valid <- RH %in% valid_RHs

    if(RH_valid) {
        result <- RH
    } else {
        cat(sprintf(
            "Invalid RH.\n  RH = %s\n",
            RH))
        result <- NA
    }
    
    return( result )
}

extract_date <- function(excel_file,
                         date_cell = "C2:C3",
                         max_date = "2020-07-01",
                         min_date = "2020-02-01"){

    date <- read_excel(excel_file,
                       range = date_cell,
                       col_types = c("date"))
    names(date) = c("date")
    date <- as_date(date$date)

    valid_date <- (date < max_date &
                   date > min_date)
    
    if(valid_date) {
        result <- date
    } else {
        cat(sprintf(
            "Invalid date. \ndate = %s \nvalid range: %s to %s \n",
            date,
            min_date,
            max_date))
        result <- NA
    }

    return( result )
}

validate_measurements <- function(measurements,
                                  valid_times){
    
    times_valid <- all(measurements$time %in% valid_times)
    temps_valid <- all(measurements$temperature_measured > 0 &
                       measurements$temperature_measured < 100)
    RHs_valid <- all(measurements$humidity_measured >= 0.0 &
                     measurements$humidity_measured <= 1.0)
    return(times_valid & temps_valid & RHs_valid)
}

extract_measurements <- function(excel_file,
                                 measurement_cells = "A60:C65",
                                 valid_times = c(-2/6, -1/12,
                                                 0, 0.5,
                                                 1, 2, 3)){
    measurements <- read_excel(excel_file,
                               range = measurement_cells)
    names(measurements) <- c("time",
                             "temperature_measured",
                             "humidity_measured")

    measurements$time[1] <- "T=-20"
    measurements$time[2] <- "T=-5"
    
    measurements$humidity_measured[1] <-
        measurements$humidity_measured[1] %>%
        str_extract("[0-9]+-") %>%
        str_extract("[0-9]+")

    ## convert humidities to decimals
    ## and times to decimal hours
    measurements <- measurements %>%
        mutate(humidity_measured = as.numeric(humidity_measured) / 100,
               time = as.numeric(substring(time,
                                           first = 3)) / 60)
    
    if(!validate_measurements(measurements,
                              valid_times)){
        cat("Invalid measurements: \n")
        print(measurements)
        result <- NA
    } else {
        result <- measurements
    }
    
    return( result )
}

extract_data <- function(excel_file){

    raw_file <- read_excel(excel_file)

    single_cell_trh <- grepl(
        "Target Chamber Humidity and Temp",
        raw_file)

    if(any(single_cell_trh)) {
        trh_col <- single_cell_trh %>%
            which() %>%
            first()
        
        trh_row <- raw_file[, trh_col] %>%
            apply(1,
                  function(x){ grepl("Target Chamber Humidity", x) }) %>%
            which() %>%
            first()

        trh_limits <- cell_limits(c(trh_row, trh_col),
                                  c(trh_row + 1, trh_col))
    } else {
        trh_limits <- NULL
    }
    
    th_row <- raw_file[, 1] %>%
        apply(1,
              function(x){ grepl("Target Humidity", x) }) %>%
        which() %>%
        first()
    
    date_row <- raw_file[,1] %>%
        apply(1,
              function(x){grepl("Date", x)}) %>%
        which() %>%
        first()

    meas_row <- raw_file[,1] %>%
        apply(1,
              function(x){grepl("AIR DRY", x)}) %>%
        which() %>%
        first()

    hum_limits <- cell_limits(c(th_row, 4), c(th_row + 1, 4))

    temp_limits <- cell_limits(c(th_row, 8), c(th_row + 1, 8))

    date_limits <- cell_limits(c(date_row, 3),
                               c(date_row + 1, 3))

    meas_limits <- cell_limits(c(meas_row + 1, 1),
                               c(meas_row + 8, 3))

    cat('extracting measurements...\n')
    measurements <- extract_measurements(
        excel_file,
        measurement_cells = meas_limits)

    cat('extracting target temperature...\n')
    T_target <- extract_target_T(excel_file,
                                 T_RH_cell = trh_limits,
                                 T_cell = temp_limits)
    
    cat('extracting target humidity...\n')
    RH_target <- extract_target_RH(excel_file,
                                   T_RH_cell = trh_limits,
                                   RH_cell = hum_limits)
    
    cat('extracting date...\n')
    date <- extract_date(excel_file,
                         date_cell = date_limits)

    cat('extracting run number...\n')
    
    run_number <- excel_file %>%
        str_extract("Run[0-9]+") %>%
        str_extract("[0-9]+") %>%
        as.numeric()
    
    result <- measurements %>%
        mutate(target_temperature = T_target,
               target_humidity = RH_target,
               date = date,
               run_number = run_number) %>%
        select(date,
               run_number,
               target_temperature,
               target_humidity,
               time,
               temperature_measured,
               humidity_measured)

    return( result )
}

main <- function(rna_metadata_path = "../dat/raw/aerosol_qpcr.csv",
                 th_metadata_directory = "../dat/raw/aerosol-metadata/",
                 outpath = "../dat/cleaned/aerosol_metadata.csv") {

    rna_meta <- read_delim(rna_metadata_path,
                           delim = ";",
                           col_types = cols()) %>%
        rename(log10_TCID50_equiv = `log10TCID50eq/ml`,
               copies_per_reaction = `copies/reaction`,
               copies_per_ml = `copies/ml`
               ) %>%
        mutate(threshold = 0.06) %>%
        select(sampling,
               CT,
               copies_per_reaction,
               copies_per_ml,
               log10_TCID50_equiv,
               run_number)
    
    directory <- th_metadata_directory
    
    excel_files <- list.files(path = directory,
                              pattern="*.xlsx",
                              full.names = TRUE)

    result <- tibble()

    for (file in excel_files){
        cat("Extracting from file ", file, "...\n")
        result <- file %>%
            extract_data() %>%
            bind_rows(result)
    }

    result <- result %>%
        distinct(run_number, .keep_all = TRUE) %>%
        arrange(run_number) %>%
        group_by(target_temperature,
                 target_humidity) %>%
        mutate(replicate = row_number()) %>%
        ungroup() %>%
        select(run_number,
               replicate) %>%
        inner_join(result, by = "run_number") %>%
        inner_join(rna_meta, by = "run_number") %>%
        arrange(run_number,
                time)
    
    write_csv(result, outpath)
}
