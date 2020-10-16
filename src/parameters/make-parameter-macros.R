#!/usr/bin/env Rscript

########################################
## filename: make-parameter-macros.R
## author: Dylan Morris <dhmorris@princeton.edu>
## create file of parameter macros
## for autopopulating priors in manuscript
#######################################

script_packages <- c(
    'stringr',    # string regex
    'dplyr'       # %>%
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
param_path <- args[1]
outpath <- args[2]

prefix <- str_match(
    outpath,
    "macros-([aA-zZ]+)-hyperparams.sty")[2] %>%
    str_to_title()

print(prefix)

sanitize <- function(param_name_string){
    return(
        param_name_string %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_remove_all(" ") %>%
        str_replace_all("[0-9]+", "Num")
    )
}

source(param_path)

macro_template <- "\\newcommand{\\%s}{%s}\n"

cat('saving macros to', outpath, '...\n')

if (file.exists(outpath))
    file.remove(outpath)

for(param_name in names(hyperparam_list)){
    cat(sprintf(macro_template,
                paste0(prefix, sanitize(param_name)),
                hyperparam_list[param_name]),
        file = outpath,
        append = TRUE)
}

warnings()
