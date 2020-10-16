  for(i_dat in 1:n_used_evap_datapoints){
    int i_class = evap_class_id[i_dat];
    real pred_mass = initial_mass[i_class] - beta[i_class] * evap_time[i_dat];
    real mu_evap = 0;
    
    if(pred_mass >= equilibrium_mass[i_class]) {
      mu_evap = pred_mass;
    } else {
      mu_evap = equilibrium_mass[i_class];
    }
    
    measured_mass[i_dat] ~ normal(mu_evap,
                                  sd_evap[i_class]);
  }


  drying_time ~ normal(mean_drying_time,
                       sd_drying_time);

  sd_evap ~ normal(mode_sd_evap,
                   sd_sd_evap);
