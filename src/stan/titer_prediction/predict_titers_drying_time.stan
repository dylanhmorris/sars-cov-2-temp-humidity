
  for(i_titer in 1:n_used_titers){
    int exp_id = titer_experiment_id[i_titer];
    int class_id = experiment_evap_class[exp_id];
    real t = titer_times[i_titer];
    real drying_t = drying_time[class_id];
    
    real ith_predicted_titer = 0;
    real ith_intercept = mean_intercept[exp_id] +
      sd_intercept[exp_id] * error_intercept[i_titer];

    if(t < drying_t) {
      ith_predicted_titer = ith_intercept - 
        transient_decay_rate[exp_id] * t;
    }
    else {
      ith_predicted_titer = ith_intercept -
        transient_decay_rate[exp_id] * drying_t -
        decay_rate[exp_id] * (t - drying_t);
    }

    // save values
    intercept[i_titer] = ith_intercept;
    sampled_titer[i_titer] = ith_predicted_titer;
}

