functions {
#include /functions/poisson_single_hit.stan
}


data {
// observed quantities
#include /general/well_data.stan

// parameters for priors  
#include /pseudoreplicates/pseudoreplicates_data.stan
#include /halflife/halflife_data.stan
  
// flags
#include /general/flags.stan

}

transformed data {}

parameters{
// this model uses non-mechanistic halflives
// with no transient phase
  
#include /halflife/halflife_parameters.stan

// this model uses data in which each point is
// an observation with a potentially different intercept
#include /pseudoreplicates/pseudoreplicates_parameters.stan

}

transformed parameters {
// declarations
#include /titer_prediction/titer_prediction_transformed_parameters.stan
#include /halflife/halflife_transformed_parameters.stan

  // calculations
#include /halflife/calculate_halflife_transformed_parameters.stan
#include /titer_prediction/predict_titers_basic.stan
}

model {

  // observation process: poisson single hit
#include /general/well_observation_process.stan
                                       

  // priors
#include /halflife/halflife_priors.stan
#include /pseudoreplicates/pseudoreplicates_priors.stan


}

generated quantities {
#include /predictive_checks/titer_predictive_checks_declarations.stan
#include /predictive_checks/titer_predictive_checks_basic.stan
}
