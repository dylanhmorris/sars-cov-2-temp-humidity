#!/usr/bin/env Rscript

########################################
## filename: table-literature-halflives.R
## author: Dylan Morris <dhmorris@princeton.edu>
## create table of half-life posterior medians
## and 95% CI for inferences from literature
#######################################

script_packages <- c(
    'virusenv',   # project functions
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'xtable',     # to create LaTeX table
    'bib2df',     # for handling meta-analysis refs
    'stringr'     # for string regex
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

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
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
lit_dat_path <- args[1]
lit_hl_path <- args[2] 
bibtex_path <- args[3] 
outpath <- args[4]

lit_dat <- read_data_for_plotting(lit_dat_path) %>%
    rename(
        DOI = study_doi)

bib_dat <- bib2df(bibtex_path)

lit_dat <- lit_dat %>%
    left_join(bib_dat,
              by = "DOI")

missing_dois <- lit_dat %>%
    filter(is.na(BIBTEXKEY))

if(dim(missing_dois)[1] > 0){
    print(missing_dois %>% distinct(DOI, .keep_all = TRUE))
    stop("Some DOIs from literature data missing from bibtex file")
}

print(lit_dat)

lit_hl_chains <- readRDS(lit_hl_path)

lit_hl_draws <- lit_hl_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(lit_dat, "experiment_id") %>%
    mutate(source = "literature")


##################################################
## calculate half-life table
##################################################

table_precisions <- c(
    0,
    0,
    0,
    0,
    0,
    0,
    2,
    2,
    2)


table_caption <- paste0(
    "\\textbf{Estimated half-lives in hours for data ",
    "from the literature, as a function of ",
    "material, temperature (T), and relative ",
    "humidity (RH).} ",
    "Estimated half-lives are reported as ",
    "posterior median and the middle 95\\% ",
    "credible interval. CCM: cell culture medium; ",
    "VTM: virus transport medium; ",
    "Resp. sec.: respiratory secretions")

table_label <- "tab:literature-halflives"

hl_table <- lit_hl_draws %>%
    mutate(hl = exp(log_half_life),
           material = ifelse(
               material == "Bulk medium",
               paste0("Bulk ", medium_group),
               material)) %>%
    filter(material != "Aerosols") %>%
    group_by(BIBTEXKEY, virus, material,
             temperature, humidity, experiment_id) %>%
    summarise(median = median(hl),
              q025 = quantile(hl, 0.025),
              q975 = quantile(hl, 0.975)) %>%
    rename(study = BIBTEXKEY) %>%
    select(-experiment_id) %>%
    mutate(study = paste0("\\citemeta{", study, "}"),
           median = as_num(median),
           q025 = as_num(q025),
           q975 = as_num(q975)) %>%
    arrange(temperature, humidity, material, study) %>%
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
## save table
####################################

cat('saving table to', outpath, '...\n')

print(hl_table,
      caption.placement = "top",
      size = "\\tiny",
      sanitize.text.function = function(x){x},
      include.rownames = FALSE,
      math.style.exponents = TRUE,
      file = outpath)

warnings()
