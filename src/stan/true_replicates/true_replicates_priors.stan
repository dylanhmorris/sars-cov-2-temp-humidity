for(i_exp in 1:n_experiments)
  true_intercept[i_exp] ~ normal(intercept_prior_mean,
                          intercept_prior_sd);
