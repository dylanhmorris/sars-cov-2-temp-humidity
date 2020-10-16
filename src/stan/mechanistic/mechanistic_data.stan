  real log_hl_20_prior_mean;
  real<lower = 0> log_hl_20_prior_sd;

  real log_k_30_20_ratio_prior_mean;
  real<lower = 0> log_k_30_20_ratio_prior_sd;
  
  real<lower = 0, upper = 1> ERH;

  real mode_sd_logconc;
  real sd_sd_logconc;

  real mean_log_shape_parameter_concavity;
  real<lower = 0> sd_log_shape_parameter_concavity;

  real mean_log_shape_parameter_steepness;
  real<lower = 0> sd_log_shape_parameter_steepness;
