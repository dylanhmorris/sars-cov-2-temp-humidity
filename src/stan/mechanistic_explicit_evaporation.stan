/*
 * filename: mechanistic_explicit_evaporation.stan
 * author: Dylan Morris <dhmorris@princetone.edu>
 * 
 * description: stan model that jointly estimates (linear
 * approximation to the) DMEM evaporation rate and the 
 * parameters of the mechanistic model of the virus 
 * inactivation rate as a function of temperature and 
 * relative humidity. This allows for inactivation during
 * the transient phase to vary over time according to the 
 * degree of evaporation that has occured.
 */


functions {
#include /functions/poisson_single_hit.stan
#include /mechanistic/mechanistic_functions.stan
}


data {
// observed quantities
#include /general/well_data.stan
#include /general/temp_humidity_data.stan
#include /evaporation_rate/evaporation_rate_data.stan
#include /mass/mass_data.stan

  
// parameters for priors
#include /pseudoreplicates/pseudoreplicates_data.stan
#include /mechanistic/mechanistic_data.stan

// flags  
#include /mechanistic/mechanistic_flags.stan  
#include /general/flags.stan
}

transformed data {
#include /evaporation_rate/evaporation_rate_transformed_data.stan
}

parameters{
#include /evaporation_rate/evaporation_rate_parameters.stan
#include /pseudoreplicates/pseudoreplicates_parameters.stan
#include /mass/mass_parameters.stan
#include /mechanistic/mechanistic_parameters.stan
}

transformed parameters {
// declare parameters

#include /titer_prediction/titer_prediction_transformed_parameters.stan
#include /mechanistic/mechanistic_transformed_parameters.stan
#include /mechanistic/mechanistic_transient_transformed_parameters.stan
#include /evaporation_rate/evaporation_rate_transformed_parameters.stan
#include /mass/mass_transformed_parameters.stan

// calculate mechanistic parameters and decay rates
#include /mechanistic/calculate_mechanistic_transformed_parameters.stan
#include /mechanistic/calculate_mechanistic_transient_transformed_parameters.stan

// predict observed titers
#include /titer_prediction/predict_titers_mechanistic_evap.stan
}

model {
  // observation processes
#include /evaporation_rate/evaporation_rate_model.stan
#include /general/well_observation_process.stan  
#include /mechanistic/concentration_observation_process.stan

  // priors
#include /mechanistic/mechanistic_priors.stan
#include /pseudoreplicates/pseudoreplicates_priors.stan
#include /mass/mass_priors.stan
}

generated quantities{
#include /predictive_checks/titer_predictive_checks_declarations.stan
#include /predictive_checks/titer_predictive_checks_evap_equation.stan
}
