if(!mass_prior_check)
  log_equilibrium_mass ~ normal(log_eq_mass_predicted,
                                sd_logconc);

if(debug){
  print("true eq mass", log_equilibrium_mass);
  print("pred eq mass", log_eq_mass_predicted);
  print("sd_logmass", sd_logconc);
  print("lposterior increment: ",
        normal_lpdf(log_equilibrium_mass | log_eq_mass_predicted,
                   sd_logconc));
  print("lposterior total: ", target());
 }
