
for(i_titer in 1:n_used_titers){
  int exp_id = titer_experiment_id[i_titer];
  int class_id = experiment_evap_class[exp_id];
  real t = titer_times[i_titer];
  
  real ith_predicted_titer = 0;
  real ith_intercept = mean_intercept[exp_id] +
    sd_intercept[exp_id] * error_intercept[i_titer];
  
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
  // note that this simplifies to
  // (initial_mass[exp_id] - equilibrium_mass[exp_id]) / beta[exp_id]
  // if we use measured concentration factors
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
  intercept[i_titer] = ith_intercept;
  sampled_titer[i_titer] = ith_predicted_titer;
 }

