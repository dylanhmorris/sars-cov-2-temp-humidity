  int<lower = 1> n_true_runs;
  int<lower = 1, upper = n_true_runs> titer_true_run_id[n_titers];

  real intercept_prior_mean;
  real intercept_prior_sd;
