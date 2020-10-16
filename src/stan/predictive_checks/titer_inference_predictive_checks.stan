  vector[n_titers] sampled_titers_pred;

  for(t_id in 1:n_titers)
    sampled_titers_pred[t_id] = normal_rng(titer_prior_mean,
                                 titer_prior_sd);
