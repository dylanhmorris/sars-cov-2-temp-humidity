#####################################
# name: Makefile
# author: Dylan Morris <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# of SARS-CoV-2 decay
# in environment
####################################

#####################################
# Directory structure
####################################

default: all

SRC := src
OUTPUT := out
DATA := dat
RAW := $(DATA)/raw
SCHEMATIC_IMAGES = $(DATA)/schematic-icons
CLEANED := $(DATA)/cleaned
MS = ms
PUBLIC = ../covid-th-public
PROJECT_PACKAGE = virusenv

# src subdirs
CLEANING_SRC = $(SRC)/cleaning
FITTING_SRC= $(SRC)/fitting
STAN_SRC = $(SRC)/stan
FIGURE_SRC = $(SRC)/figures
TABLE_SRC = $(SRC)/tables
PARAMS := $(SRC)/parameters

# output and manuscript subdirs
MCMC_CHAINS := $(OUTPUT)/mcmc_chains
FIGURE_DIR = $(OUTPUT)/figures
TABLE_DIR = $(OUTPUT)/tables

PRIOR_CHECK_DIR_NAME = prior-checks
POSTERIOR_CHECK_DIR_NAME = posterior-checks

MS_FIG_DIR = $(MS)/figures
MS_TABLE_DIR = $(MS)/tables
MS_BIBTEX_FILE = $(MS)/bibliography.bib
MS_MACRO_DIR = $(MS)/macros

#####################################
# File extensions and the like
####################################
CHAINS_SUFFIX = -chains.Rds
PRIOR_CHECK_SUFFIX = -prior-check.Rds

#####################################
# Expected bash settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

MKDIR := @mkdir -p
CP := @cp
RM := rm -rf

R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)

# R creates a blank Rplots.pdf when run
# from the command line. This removes it.
FIG_CLEANUP = @$(RM) Rplots.pdf

#####################################
# Installation / dependencies
#
# Rules for prepping analysis
#####################################

.PHONY: depend

depend:
	$(R_COMMAND) $(SRC)/install_needed_packages.R


################################
# model names
################################

# half-life models / generic names
PLA_NAME = plastic
EVAP_NAME = evaporation
EVAP_SALT_NAME = evaporation-modeled

LIT_NAME = literature

SARS_MERS_NAME = sars-mers

HL_MODEL_NAMES = $(PLA_NAME) $(SARS_MERS_NAME) 

# mechanistic kinetics models
MECH_$(EVAP_NAME)_NAME = mechanistic-$(EVAP_NAME)
MECH_$(EVAP_SALT_NAME)_NAME = mechanistic-$(EVAP_SALT_NAME)

MECH_MODEL_NAMES = $(MECH_$(EVAP_NAME)_NAME) $(MECH_$(EVAP_SALT_NAME)_NAME)
# titer inference models
INFER_$(PLA_NAME)_NAME = inferred-$(PLA_NAME)-titers
INFER_$(SARS_MERS_NAME)_NAME = inferred-$(SARS_MERS_NAME)-titers

TITER_MODEL_NAMES = $(INFER_$(PLA_NAME)_NAME) $(INFER_$(SARS_MERS_NAME)_NAME)

####################
# groups of models
####################
# models with a evaporation phase
EVAPORATION_MODEL_NAMES =  $(HL_MODEL_NAMES) $(MECH_MODEL_NAMES)
# all models except raw titer inference
MODEL_NAMES = $(HL_MODEL_NAMES) $(LIT_NAME) $(MECH_MODEL_NAMES)
# all models
MODEL_NAMES_INCL_TITER = $(MODEL_NAMES) $(TITER_MODEL_NAMES)

# titers for additional models
INFER_$(MECH_$(PLA_NAME)_NAME)_NAME = inferred-$(PLA_NAME)-titers
INFER_$(MECH_$(SALT_NAME)_NAME)_NAME = inferred-$(PLA_NAME)-titers
INFER_$(MECH_$(EVAP_NAME)_NAME)_NAME = inferred-$(PLA_NAME)-titers
INFER_$(MECH_$(EVAP_SALT_NAME)_NAME)_NAME = inferred-$(PLA_NAME)-titers


#####################################
# data locations
#####################################
$(PLA_NAME)_DATAFILE = $(PLA_NAME)-data.csv
$(LIT_NAME)_DATAFILE = literature-data.csv
$(SARS_MERS_NAME)_DATAFILE = $(SARS_MERS_NAME)-data.csv
EVAPORATION_DATAFILE = dmem-evaporation-data.csv
EVAP_ENDPOINTS_FILE = dmem-evap-endpoints.csv

SCHEMATIC_FILENAMES = effloresced-droplet.png evaporating-droplet.png large-droplet.png small-droplet.png

SCHEMATIC_PATHS = $(addprefix $(SCHEMATIC_IMAGES)/, $(SCHEMATIC_FILENAMES))


RAW_TH_DATA = $(RAW)/$(TH_DATAFILE)

## shortnames for cleaned data
$(PLA_NAME)_DATA = $(CLEANED)/$($(PLA_NAME)_DATAFILE)
$(LIT_NAME)_DATA = $(CLEANED)/$($(LIT_NAME)_DATAFILE)
$(SARS_MERS_NAME)_DATA = $(CLEANED)/$($(SARS_MERS_NAME)_DATAFILE)
EVAPORATION_DATA = $(CLEANED)/$(EVAPORATION_DATAFILE)

CLEANED_DATA = $($(PLA_NAME)_DATA) $($(LIT_NAME)_DATA) $($(SARS_MERS_NAME)_DATA) $(EVAPORATION_DATA)


###########################
## associate data to models

$(PLA_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA) $(EVAPORATION_DATA)
$(LIT_NAME)_FITTING_DATA = $($(LIT_NAME)_DATA)
$(SARS_MERS_NAME)_FITTING_DATA = $($(SARS_MERS_NAME)_DATA) $(EVAPORATION_DATA)
$(EVAP_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA) $(EVAPORATION_DATA)
$(EVAP_SALT_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA) $(EVAPORATION_DATA)


$(MECH_$(PLA_NAME)_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA)
$(MECH_$(SALT_NAME)_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA)

$(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA) $(EVAPORATION_DATA)
$(MECH_$(EVAP_SALT_NAME)_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA) $(EVAPORATION_DATA)

$(INFER_$(PLA_NAME)_NAME)_FITTING_DATA = $($(PLA_NAME)_DATA)
$(INFER_$(SARS_MERS_NAME)_NAME)_FITTING_DATA = $($(SARS_MERS_NAME)_DATA)
#####################################
# Code locations
#####################################

## enumerate data cleaning scripts
TH_CLEANING_SCRIPT = $(CLEANING_SRC)/clean-th-data.R
LIT_CLEANING_SCRIPT = $(CLEANING_SRC)/clean-literature-data.R


## enumerate model fitting scripts
PLA_FITTING_SCRIPT = $(FITTING_SRC)/fit_plastic_evap_model.R
INFER_FITTING_SCRIPT = $(FITTING_SRC)/fit_plastic_model.R
LIT_FITTING_SCRIPT = $(FITTING_SRC)/fit_literature_model.R

$(EVAP_NAME)_FITTING_SCRIPT = $(FITTING_SRC)/fit_plastic_evap_model.R

## associate models to fitting scripts
$(PLA_NAME)_FITTING_SCRIPT = $(PLA_FITTING_SCRIPT)
$(LIT_NAME)_FITTING_SCRIPT = $(LIT_FITTING_SCRIPT)
$(SARS_MERS_NAME)_FITTING_SCRIPT = $(PLA_FITTING_SCRIPT)

$(MECH_$(EVAP_NAME)_NAME)_FITTING_SCRIPT = $($(EVAP_NAME)_FITTING_SCRIPT)

$(INFER_$(PLA_NAME)_NAME)_FITTING_SCRIPT = $(INFER_FITTING_SCRIPT)
$(INFER_$(SARS_MERS_NAME)_NAME)_FITTING_SCRIPT = $(INFER_FITTING_SCRIPT)

## stan source files
$(PLA_NAME)_MODEL_SRC = halflives_explicit_evaporation.stan
$(LIT_NAME)_MODEL_SRC = halflives_fraction_viable.stan
$(SARS_MERS_NAME)_MODEL_SRC = halflives_explicit_evaporation.stan
$(EVAP_NAME)_MODEL_SRC = evaporation_rate.stan
$(MECH_$(EVAP_NAME)_NAME)_MODEL_SRC = mechanistic_explicit_evaporation.stan

$(INFER_$(PLA_NAME)_NAME)_MODEL_SRC = infer_titers.stan
$(INFER_$(SARS_MERS_NAME)_NAME)_MODEL_SRC = infer_titers.stan

## chain diagnostics scripts
DIAGNOSTIC_SCRIPT = $(SRC)/chain_diagnostics.R

#####################################
# parameter locations
#####################################

## paths to hyperparameters
PLA_MECH_HYPERS = $(PARAMS)/plastic_hyperparams.R
LIT_HYPERS = $(PARAMS)/literature_hyperparams.R
EVAP_HYPERS = $(PARAMS)/plastic_hyperparams.R
INFER_HYPERS = $(PARAMS)/infer_hyperparams.R

$(EVAP_NAME)_HYPERS = $(EVAP_HYPERS)
$(PLA_NAME)_HYPERS = $(PLA_MECH_HYPERS)
$(LIT_NAME)_HYPERS = $(LIT_HYPERS)
$(SARS_MERS_NAME)_HYPERS = $(PLA_MECH_HYPERS)

$(MECH_$(PLA_NAME)_NAME)_HYPERS = $(PLA_MECH_HYPERS)
$(MECH_$(SALT_NAME)_NAME)_HYPERS = $(PLA_MECH_HYPERS)
$(MECH_$(EVAP_NAME)_NAME)_HYPERS = $(EVAP_HYPERS)

$(INFER_$(PLA_NAME)_NAME)_HYPERS = $(INFER_HYPERS)
$(INFER_$(SARS_MERS_NAME)_NAME)_HYPERS = $(INFER_HYPERS)

#####################################
# mcmc output locations
#####################################

## half lives
HL_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(addsuffix $(CHAINS_SUFFIX), $(HL_MODEL_NAMES)))

## mechanistic kinetics
MECH_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(addsuffix $(CHAINS_SUFFIX), $(MECH_MODEL_NAMES)))

## literature meta-analysis
LIT_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(addsuffix $(CHAINS_SUFFIX), $(LIT_NAME)))

## titer inference
TITER_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(addsuffix $(CHAINS_SUFFIX), $(TITER_MODEL_NAMES)))

## evaporation inference
EVAP_CHAINS = $(MCMC_CHAINS)/$(EVAP_NAME)$(CHAINS_SUFFIX)

#############################
## prior predictive checks
#############################

PRIOR_CHECK_NAMES = $(addsuffix $(PRIOR_CHECK_SUFFIX), $(MODEL_NAMES_INCL_TITER))
PRIOR_CHECK_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(PRIOR_CHECK_NAMES))


#############################
## all mcmc chains
#############################

CHAIN_PATHS = $(EVAP_CHAINS) $(PRIOR_CHECK_CHAINS) $(HL_CHAINS) $(MECH_CHAINS) $(LIT_CHAINS) $(TITER_CHAINS) 

CHAIN_DIAGNOSTICS = $(OUTPUT)/chain_diagnostics.csv


#############################
## all figures
#############################

PRIOR_CHECK_FIGURES = $(addprefix $(PRIOR_CHECK_DIR_NAME)/figure-, $(addsuffix -prior-check.pdf, $(MODEL_NAMES_INCL_TITER))) $(addprefix $(PRIOR_CHECK_DIR_NAME)/figure-evap-phase-, $(addsuffix -prior-check.pdf, $(EVAPORATION_MODEL_NAMES)))

POSTERIOR_CHECK_FIGURES = $(addprefix $(POSTERIOR_CHECK_DIR_NAME)/figure-, $(addsuffix -posterior-check.pdf, $(MODEL_NAMES))) $(addprefix $(POSTERIOR_CHECK_DIR_NAME)/figure-evap-phase-, $(addsuffix -posterior-check.pdf, $(EVAPORATION_MODEL_NAMES)))

FIGURES = figure-mech-fit.pdf figure-mech-fit-evap-phase.pdf figure-mech-fit-mod-conc.pdf figure-mech-fit-mod-conc-evap-phase.pdf figure-halflives-schematic.pdf figure-extrapolation.pdf figure-concentration-factor-mech.pdf figure-concentration-factor-model.pdf figure-predict-literature.pdf figure-hl-fit.pdf figure-hl-fit-evap-phase.pdf figure-evaporation.pdf figure-prisma.pdf figure-sars-mers-regression.pdf $(PRIOR_CHECK_FIGURES) $(POSTERIOR_CHECK_FIGURES)

FIGURE_PATHS = $(addprefix $(FIGURE_DIR)/, $(FIGURES))
MS_FIGURES = $(addprefix $(MS_FIG_DIR)/, $(FIGURES))

#############################
## all tables
#############################

TABLES = table-$(PLA_NAME)-halflives.tex table-$(SARS_MERS_NAME)-halflives.tex table-$(LIT_NAME)-halflives.tex

TABLE_PATHS = $(addprefix $(TABLE_DIR)/, $(TABLES))
MS_TABLES = $(addprefix $(MS_TABLE_DIR)/, $(TABLES))

#############################
## all macros
#############################

MACRO_NAMES = macros-$(PLA_NAME)-halflives.sty macros-$(SARS_MERS_NAME)-halflives.sty macros-plastic-hyperparams.sty macros-literature-hyperparams.sty macros-infer-hyperparams.sty
MACRO_PATHS = $(addprefix $(MS_MACRO_DIR)/, $(MACRO_NAMES))


#####################################
# Rules
#
# definition of dependency
# tree and specification of
# rules for doing stuff
#####################################

##########################
# rules for data cleaning
##########################

$(CLEANED)/$(EVAPORATION_DATAFILE): $(CLEANING_SRC)/clean-dmem-evap-data.R $(RAW)/$(EVAPORATION_DATAFILE) $(RAW)/$(EVAP_ENDPOINTS_FILE)
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@

$(CLEANED)/literature-data.csv: $(CLEANING_SRC)/clean-literature-data.R $(RAW)/literature-data.csv
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@

$(CLEANED)/sars-mers-data.csv: $(CLEANING_SRC)/clean-plastic-data.R $(RAW)/th-data.csv
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@

$(CLEANED)/%-data.csv: $(CLEANING_SRC)/clean-%-data.R $(RAW)/th-data.csv
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@


#####################################
# rules for model fitting and post-processing
#####################################

.SECONDEXPANSION:

$(MCMC_CHAINS)/$(MECH_$(EVAP_SALT_NAME)_NAME)$(CHAINS_SUFFIX): $($(MECH_$(EVAP_NAME)_NAME)_FITTING_SCRIPT) $(STAN_SRC)/$($(MECH_$(EVAP_NAME)_NAME)_MODEL_SRC) $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $($(MECH_$(EVAP_NAME)_NAME)_HYPERS)
	@echo Making evap salt: $@... "\n\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ FALSE FALSE

$(MCMC_CHAINS)/$(MECH_$(EVAP_SALT_NAME)_NAME)$(PRIOR_CHECK_SUFFIX): $($(MECH_$(EVAP_NAME)_NAME)_FITTING_SCRIPT) $(STAN_SRC)/$($(MECH_$(EVAP_NAME)_NAME)_MODEL_SRC) $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $($(MECH_$(EVAP_NAME)_NAME)_HYPERS)
	@echo Making evap salt prior check: $@... "\n\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ TRUE FALSE


# generic model fitting rule
$(MCMC_CHAINS)/%$(CHAINS_SUFFIX): $$($$*_FITTING_SCRIPT) $(STAN_SRC)/$$($$*_MODEL_SRC) $$($$*_FITTING_DATA) $$($$*_HYPERS)
	@echo Making $@... "\n\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ FALSE TRUE

$(MCMC_CHAINS)/%$(PRIOR_CHECK_SUFFIX): $$($$*_FITTING_SCRIPT) $(STAN_SRC)/$$($$*_MODEL_SRC) $$($$*_FITTING_DATA) $$($$*_HYPERS)
	@echo "\nMaking" $@... "\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ TRUE TRUE

###########################
## diagnostics for MCMC
###########################

$(CHAIN_DIAGNOSTICS): $(DIAGNOSTIC_SCRIPT) $(CHAIN_PATHS)
	$(MKDIR) $(OUTPUT)
	$(R_COMMAND) $^ $@



#####################################
# rules for figure generation
#####################################

## useful shorthand
INF_PLA_CHAINS = $(MCMC_CHAINS)/$(INFER_$(PLA_NAME)_NAME)$(CHAINS_SUFFIX)
MECH_EVAP_CHAINS = $(MCMC_CHAINS)/$(MECH_$(EVAP_NAME)_NAME)$(CHAINS_SUFFIX)
MECH_EVAP_SALT_CHAINS = $(MCMC_CHAINS)/$(MECH_$(EVAP_SALT_NAME)_NAME)$(CHAINS_SUFFIX)

$(FIGURE_DIR)/figure-halflives-schematic.pdf: $(FIGURE_SRC)/figure-halflives-schematic.R $($(PLA_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(INFER_$(PLA_NAME)_NAME)$(CHAINS_SUFFIX) $(SCHEMATIC_PATHS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-mech-fit.pdf: $(FIGURE_SRC)/figure-regression.R $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $(MECH_EVAP_CHAINS) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-extrapolation.pdf: $(FIGURE_SRC)/figure-extrapolation.R $($(PLA_NAME)_FITTING_DATA) $($(LIT_NAME)_FITTING_DATA) $($(SARS_MERS_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(LIT_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(SARS_MERS_NAME)$(CHAINS_SUFFIX) $(MECH_EVAP_CHAINS) $(MECH_EVAP_SALT_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-extrapolation.jpg: $(FIGURE_SRC)/figure-extrapolation.R $($(PLA_NAME)_FITTING_DATA) $($(LIT_NAME)_FITTING_DATA) $($(SARS_MERS_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(LIT_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(SARS_MERS_NAME)$(CHAINS_SUFFIX) $(MECH_EVAP_CHAINS) $(MECH_EVAP_SALT_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(FIGURE_DIR)/figure-mech-fit-evap-phase.pdf: $(FIGURE_SRC)/figure-regression.R $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $(MECH_EVAP_CHAINS) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-mech-fit-mod-conc.pdf: $(FIGURE_SRC)/figure-regression.R $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $(MECH_EVAP_SALT_CHAINS) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-mech-fit-mod-conc-evap-phase.pdf: $(FIGURE_SRC)/figure-regression.R $($(MECH_$(EVAP_NAME)_NAME)_FITTING_DATA) $(MECH_EVAP_SALT_CHAINS) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(FIGURE_DIR)/figure-hl-fit.pdf: $(FIGURE_SRC)/figure-regression.R $($(PLA_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-hl-fit-evap-phase.pdf: $(FIGURE_SRC)/figure-regression.R $($(PLA_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(INF_PLA_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-sars-mers-regression.pdf: $(FIGURE_SRC)/figure-sars-mers-regression.R $($(SARS_MERS_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(SARS_MERS_NAME)$(CHAINS_SUFFIX)  $(MCMC_CHAINS)/$(INFER_$(SARS_MERS_NAME)_NAME)$(CHAINS_SUFFIX)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-predict-literature.pdf: $(FIGURE_SRC)/figure-predict-literature.R $($(PLA_NAME)_FITTING_DATA) $($(LIT_NAME)_FITTING_DATA) $($(SARS_MERS_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(LIT_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(SARS_MERS_NAME)$(CHAINS_SUFFIX) $(MECH_EVAP_SALT_CHAINS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-concentration-factor-mech.pdf: $(FIGURE_SRC)/figure-concentration-factor.R $($(PLA_NAME)_FITTING_DATA) $(MECH_EVAP_CHAINS) 
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-concentration-factor-model.pdf: $(FIGURE_SRC)/figure-concentration-factor.R $($(PLA_NAME)_FITTING_DATA) $(MECH_EVAP_SALT_CHAINS) 
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-concentration-factor-mech-check.pdf: $(FIGURE_SRC)/figure-concentration-factor.R $($(PLA_NAME)_DATA) $(MCMC_CHAINS)/$(MECH_$(EVAP_NAME)_NAME)$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-concentration-factor-mod-check.pdf: $(FIGURE_SRC)/figure-concentration-factor.R $($(PLA_NAME)_DATA) $(MCMC_CHAINS)/$(MECH_$(EVAP_SALT_NAME)_NAME)$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure-evaporation.pdf: $(FIGURE_SRC)/figure-evaporation.R $($(PLA_NAME)_DATA) $(EVAPORATION_DATA) $(MCMC_CHAINS)/$(PLA_NAME)$(CHAINS_SUFFIX)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(PRIOR_CHECK_DIR_NAME)/figure-literature-prior-check.pdf: $(FIGURE_SRC)/figure-literature-predictive-check.R $($(LIT_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(LIT_NAME)$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(POSTERIOR_CHECK_DIR_NAME)/figure-literature-posterior-check.pdf: $(FIGURE_SRC)/figure-literature-predictive-check.R $($(LIT_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(LIT_NAME)$(CHAINS_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(PRIOR_CHECK_DIR_NAME)/figure-evap-phase-%-prior-check.pdf: $(FIGURE_SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$(EVAP_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/%$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(PRIOR_CHECK_DIR_NAME)/figure-inferred-%-titers-prior-check.pdf: $(FIGURE_SRC)/figure-inference-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(FIGURE_DIR)/$(PRIOR_CHECK_DIR_NAME)/figure-%-prior-check.pdf: $(FIGURE_SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$(EVAP_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/%$(PRIOR_CHECK_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(POSTERIOR_CHECK_DIR_NAME)/figure-evap-phase-%-posterior-check.pdf: $(FIGURE_SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/%$(CHAINS_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/$(POSTERIOR_CHECK_DIR_NAME)/figure-%-posterior-check.pdf: $(FIGURE_SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$$(INFER_$$*_NAME)$(CHAINS_SUFFIX) $(MCMC_CHAINS)/%$(CHAINS_SUFFIX)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

## rule to copy figures from output directory
## to the manuscript directory
$(FIGURE_DIR)/figure-prisma.pdf: $(DATA)/figure-prisma.pdf
	$(MKDIR) $(dir $@)
	$(CP) $< $@

$(FIGURE_DIR)/%.pdf: $(FIGURE_SRC)/%.R
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

## rule to copy figures from output directory
## to the manuscript directory
$(MS_FIG_DIR)/%.pdf: $(FIGURE_DIR)/%.pdf
	$(MKDIR) $(dir $@)
	$(CP) $< $@


#####################################
# rules for tables
#####################################


$(TABLE_DIR)/table-literature-halflives.tex: $(TABLE_SRC)/table-literature-halflives.R $($(LIT_NAME)_FITTING_DATA) $(MCMC_CHAINS)/$(LIT_NAME)$(CHAINS_SUFFIX) $(MS_BIBTEX_FILE)
	$(MKDIR) $(TABLE_DIR)
	$(R_COMMAND) $^ $@

$(TABLE_DIR)/table-%-halflives.tex: $(TABLE_SRC)/table-halflives.R  $$($$*_FITTING_DATA) $(MCMC_CHAINS)/%$(CHAINS_SUFFIX)
	$(MKDIR) $(TABLE_DIR)
	$(R_COMMAND) $^ $@

$(TABLE_DIR)/%.tex: $(TABLE_SRC)/%.R
	$(MKDIR) $(TABLE_DIR)
	$(R_COMMAND) $^ $@

## rule to copy tables from output directory
## to the manuscript directory
$(MS_TABLE_DIR)/%.tex: $(TABLE_DIR)/%.tex
	$(MKDIR) $(MS_TABLE_DIR)
	$(CP) $< $@

#####################################
# rules for dynamic macros
#####################################

$(MS_MACRO_DIR)/macros-%-halflives.sty: $(TABLE_SRC)/table-halflives.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/%$(CHAINS_SUFFIX)
	$(MKDIR) $(MS_MACRO_DIR)
	$(R_COMMAND) $^ $@

$(MS_MACRO_DIR)/macros-%-hyperparams.sty: $(PARAMS)/make-parameter-macros.R $(PARAMS)/%_hyperparams.R
	$(MKDIR) $(MS_MACRO_DIR)
	$(R_COMMAND) $^ $@

#####################################
# convenience rules for making
# various quantities
#####################################
.PHONY: data
data: $(CLEANED_DATA)

.PHONY: chains titer_chains hl_chains mech_chains prior_check
chains: $(CHAIN_PATHS)
hl_chains: $(HL_CHAINS)
mech_chains: $(MECH_CHAINS)
titer_chains: $(TITER_CHAINS) 
prior_checks: $(PRIOR_CHECK_CHAINS)

.PHONY: diagnostics
diagnostics: $(CHAIN_DIAGNOSTICS)

.PHONY: figures
figures: $(FIGURE_PATHS) $(MS_FIGURES)

.PHONY: tables
tables: $(TABLE_PATHS) $(MS_TABLES)

.PHONY: macros
macros: $(MACRO_PATHS)

## get info about paths
.PHONY: echo_models echo_figures echo_chains echo_checks

echo_models:
	@echo $(MODEL_NAMES_INCL_TITER)
	@echo $(CHAIN_PATHS)

echo_figures:
	@echo $(FIGURE_PATHS)

echo_chains:
	@echo $(CHAIN_PATHS)

echo_checks:
	@echo $(PRIOR_CHECK_CHAINS)

echo_check_figs:
	@echo $(PRIOR_CHECK_FIGURES)


## cleanup rules
.PHONY: deltemp delcached clean

## remove emacs tempfiles etc.
deltemp:
	$(RM) src/*~*
	$(RM) src/*#*
	$(RM) src/*/*~*
	$(RM) src/*/*#*
	$(RM) src/stan/*/*~*
	$(RM) src/stan/*/*#*

## remove cached stan models
delcached:
	$(RM) src/stan/*.rds
	$(RM) src/stan/*.Rds

## full clean
clean: deltemp delcached
	$(RM) $(OUTPUT)
	$(RM) $(CLEANED)
	$(RM) $(PUBLIC)

## produce public copy
.PHONY: public

public: depend deltemp
	$(RM) $(PUBLIC)
	$(MKDIR) $(PUBLIC)
	$(MKDIR) $(PUBLIC)/$(OUTPUT)
	$(CP) -r $(DATA) $(PUBLIC)/$(DATA)
	$(CP) -r $(FIGURE_DIR) $(PUBLIC)/$(FIGURE_DIR)
	$(CP) -r $(TABLE_DIR) $(PUBLIC)/$(TABLE_DIR)
	$(CP) -r $(SRC) $(PUBLIC)/$(SRC)
	$(CP) -r $(PROJECT_PACKAGE) $(PUBLIC)/$(PROJECT_PACKAGE)
	$(CP) .gitignore $(PUBLIC)/.gitignore
	$(CP) README.md $(PUBLIC)/README.md
	$(CP) LICENSE.txt $(PUBLIC)/LICENSE.txt
	$(CP) Makefile $(PUBLIC)/Makefile

## make everything
all: depend data chains diagnostics figures tables macros

