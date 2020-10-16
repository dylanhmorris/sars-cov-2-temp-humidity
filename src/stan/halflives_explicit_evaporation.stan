functions {
#include /functions/poisson_single_hit.stan
}

data {
// observed quantities
#include /general/well_data.stan

// parameters for priors  
#include /halflife/halflife_data.stan
#include /halflife/halflife_transient_data.stan
#include /pseudoreplicates/pseudoreplicates_data.stan
#include /evaporation_rate/evaporation_rate_data.stan

// flags
#include /general/flags.stan
  
}

transformed data {
#include /evaporation_rate/evaporation_rate_transformed_data.stan
}



parameters{
  // this model explicitly estimates evaporation rates
  // and drying times

#include /evaporation_rate/evaporation_rate_parameters.stan

// this model uses non-mechanistic halflives
// with a transient phase
#include /halflife/halflife_parameters.stan
#include /halflife/halflife_transient_parameters.stan
  
// this model uses data in which each point is
// an observation with a potentially different intercept
#include /pseudoreplicates/pseudoreplicates_parameters.stan

}

transformed parameters {
  // declarations
#include /titer_prediction/titer_prediction_transformed_parameters.stan
#include /halflife/halflife_transformed_parameters.stan
#include /halflife/halflife_transient_transformed_parameters.stan
#include /evaporation_rate/evaporation_rate_transformed_parameters.stan

  // calculations
#include /halflife/calculate_halflife_transformed_parameters.stan
#include /halflife/calculate_halflife_transient_transformed_parameters.stan

#include /titer_prediction/predict_titers_drying_time.stan
}

model {
  // observation processes
#include /evaporation_rate/evaporation_rate_model.stan
#include /general/well_observation_process.stan
                                       
  // priors
#include /halflife/halflife_priors.stan
#include /halflife/halflife_transient_priors.stan
#include /pseudoreplicates/pseudoreplicates_priors.stan
}

generated quantities {
#include /predictive_checks/titer_predictive_checks_declarations.stan
#include /predictive_checks/titer_predictive_checks_drying.stan
}
