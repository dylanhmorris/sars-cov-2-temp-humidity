  mean_intercept ~ normal(intercept_prior_mean,
                          intercept_prior_sd);
  

  error_intercept ~ normal(0, 1);
  
  sd_intercept  ~ normal(mode_sd_intercept,
                         sd_sd_intercept);
