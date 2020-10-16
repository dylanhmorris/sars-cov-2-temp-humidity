  log_hl_20_dry ~ normal(log_hl_20_prior_mean,
                         log_hl_20_prior_sd);
  log_hl_20_wet ~ normal(log_hl_20_prior_mean,
                         log_hl_20_prior_sd);

  log_k_30_20_ratio_dry ~ normal(log_k_30_20_ratio_prior_mean,
                                 log_k_30_20_ratio_prior_sd);
  log_k_30_20_ratio_wet ~ normal(log_k_30_20_ratio_prior_mean,
                                 log_k_30_20_ratio_prior_sd);

  log_shape_parameter_concavity ~ normal(mean_log_shape_parameter_concavity,
                                         sd_log_shape_parameter_concavity);
  log_shape_parameter_steepness ~ normal(mean_log_shape_parameter_steepness,
                                         sd_log_shape_parameter_steepness);
  sd_logconc ~ normal(mode_sd_logconc,
                      sd_sd_logconc);
