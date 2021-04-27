functions {
#include /functions/poisson_single_hit.stan
}

data {
// observed quantities
#include /general/well_data.stan

// parameters for priors  
#include /true_replicates/true_replicates_data.stan
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
#include /true_replicates/true_replicates_parameters.stan

}

transformed parameters {
// declarations
#include /titer_prediction/titer_prediction_transformed_parameters.stan
#include /halflife/halflife_transformed_parameters.stan

  // calculations
#include /halflife/calculate_halflife_transformed_parameters.stan
#include /titer_prediction/predict_titers_true_replicate.stan
}

model {

  // observation process: poisson single hit
#include /general/well_observation_process.stan
                                       

  // priors
#include /halflife/halflife_priors.stan
#include /true_replicates/true_replicates_priors.stan


}

generated quantities {
}
