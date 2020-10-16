functions {}

data {
  int<lower = 1> n_experiments;
#include /evaporation_rate/evaporation_rate_data.stan    
}

transformed data {
#include /evaporation_rate/evaporation_rate_transformed_data.stan    
}

parameters{
#include /evaporation_rate/evaporation_rate_parameters.stan  
}

transformed parameters {
#include /evaporation_rate/evaporation_rate_transformed_parameters.stan  
}

model {
#include /evaporation_rate/evaporation_rate_model.stan
}

generated quantities{}
