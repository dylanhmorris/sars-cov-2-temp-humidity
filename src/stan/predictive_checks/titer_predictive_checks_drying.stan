
  for(i_titer in 1:n_titers){
    int i_exp = titer_experiment_id[i_titer];
    int class_id = experiment_evap_class[i_exp];

    real t = titer_times[i_titer];

    real drying_t = drying_time[class_id];
    
    real ith_predicted_titer = 0;
    real ith_intercept = normal_rng(mean_intercept[i_exp],
                                    sd_intercept[i_exp]);
    
    if(t < drying_t) {
      ith_predicted_titer = ith_intercept - 
        transient_decay_rate[i_exp] * t;
    }
    else {
      ith_predicted_titer = ith_intercept -
        transient_decay_rate[i_exp] * drying_t -
        decay_rate[i_exp] * (t - drying_t);
    }      
    
    // save values
    intercepts_pred[i_titer] = ith_intercept;
    sampled_titers_pred[i_titer] = ith_predicted_titer;
  }
