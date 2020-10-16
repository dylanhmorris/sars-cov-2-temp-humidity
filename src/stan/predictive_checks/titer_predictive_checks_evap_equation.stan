
for(i_titer in 1:n_titers){
  int exp_id = titer_experiment_id[i_titer];
  real t = titer_times[i_titer];
  int class_id = experiment_evap_class[exp_id];

  real ith_predicted_titer = 0;
  real ith_intercept = normal_rng(mean_intercept[exp_id],
                                  sd_intercept[exp_id]);
  
  real eq_conc = 0;
  real drying_t = 0;
  real B = beta[class_id] / (initial_mass[class_id] -
                           solute_mass_init[class_id]);
  
  if(empirical_salt){ 
    eq_conc =
         exp(log_concentration_factor[class_id]);
  } else {
    eq_conc =
      exp(log_concentration_factor_predicted[class_id]);
  }
  drying_t = (1 - (1 / eq_conc)) / B;
    
  if(t < drying_t) {
    ith_predicted_titer = ith_intercept +
      (transient_decay_rate[exp_id] / B) *
      log10(1 - B * t);
  } else {
    ith_predicted_titer = ith_intercept +
      (transient_decay_rate[exp_id] / B) *
      log10(1 - B * drying_t) -
      decay_rate[exp_id] * (t - drying_t);
  }

    // save values
  intercepts_pred[i_titer] = ith_intercept;
  sampled_titers_pred[i_titer] = ith_predicted_titer;
 }

