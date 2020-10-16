  int<lower = 1> n_transient_classes;
  int<lower = 1, upper = n_transient_classes> transient_class_id[n_experiments];
  real log_hl_transient_prior_mean;
  real<lower = 0> log_hl_transient_prior_sd;
