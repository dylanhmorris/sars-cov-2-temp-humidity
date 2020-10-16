  vector[n_experiments] log_concentration_factor_predicted;
  vector[n_experiments] log_eq_mass_predicted;

  vector[n_experiments] decay_rate;  

  real<lower = 0> k_20_dry;
  real<lower = 0> k_20_wet;
  
  real log_A_dry;
  real log_A_wet;

  real<lower = 0> E_a_dry;
  real<lower = 0> E_a_wet;
