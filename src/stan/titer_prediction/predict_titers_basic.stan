  for(i_titer in 1:n_used_titers){
    int exp_id = titer_experiment_id[i_titer];
    real t = titer_times[i_titer];
    
    real ith_predicted_titer = 0;
    real ith_intercept = mean_intercept[exp_id] +
      sd_intercept[exp_id] * error_intercept[i_titer];
    
    ith_predicted_titer = ith_intercept - 
      decay_rate[exp_id] * t;

    // save values
    intercept[i_titer] = ith_intercept;
    sampled_titer[i_titer] = ith_predicted_titer;
}
