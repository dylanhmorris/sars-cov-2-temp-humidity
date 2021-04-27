  /* filename: infer_titers.stan
   * author: Dylan H. Morris <dylan@dylanhmorris.com>
   * 
   * description: estimate virus titers 
   * (log10 TCID50 / vol) directly 
   * from raw titration (well) data
   */


functions {
#include /functions/poisson_single_hit.stan
}


data {
#include /general/well_data.stan
#include /general/flags.stan

  // priors
#include /titer_inference/titer_inference_data.stan
}

transformed data {}


parameters{
#include /titer_inference/titer_inference_parameters.stan
}

transformed parameters {}

model {

  // observation process
#include /general/well_observation_process.stan  


  // priors
#include /titer_inference/titer_inference_priors.stan
}

generated quantities {
#include /predictive_checks/titer_inference_predictive_checks.stan
}
