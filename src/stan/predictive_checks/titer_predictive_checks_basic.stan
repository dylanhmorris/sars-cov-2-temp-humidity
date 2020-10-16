
  for(i_titer in 1:n_titers){
    int i_exp = titer_experiment_id[i_titer];
    real t = titer_times[i_titer];
    
    real ith_predicted_titer = 0;
    real ith_intercept = normal_rng(mean_intercept[i_exp],
                                    sd_intercept[i_exp]);
    
    ith_predicted_titer = ith_intercept - 
      decay_rate[i_exp] * t;
      
    // save values
    intercepts_pred[i_titer] = ith_intercept;
    sampled_titers_pred[i_titer] = ith_predicted_titer;
  }
