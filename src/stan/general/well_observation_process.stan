  for (i_dat in 1:n_used_datapoints) {
    int i_titer = titer_id[i_dat];
    real dilute_dose = sampled_titer[i_titer] + dilution[i_dat];
    real virions = log(2) * pow(10, dilute_dose);

    well_status[i_dat] ~ poisson_single_hit(virions);
  }
