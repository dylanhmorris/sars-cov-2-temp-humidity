// calculate key mechanistic parameters
  E_a_dry = R() * log_k_30_20_ratio_dry / temp_inv_diff(20, 30);
  E_a_wet = R() * log_k_30_20_ratio_wet / temp_inv_diff(20, 30);
  
  k_20_dry = log10(2) / exp(log_hl_20_dry);
  k_20_wet = log10(2) / exp(log_hl_20_wet);
  
  
  log_A_dry = log_arrenhius_A(k_20_dry,
                              E_a_dry,
                              to_Kelvin(20));


  log_A_wet = log_arrenhius_A(k_20_wet,
                              E_a_wet,
                              to_Kelvin(20));

  // calculate decay rates for each experiment 
  for(i_exp in 1:n_experiments) {
    real RH = rel_humidity[i_exp];
    real T = to_Kelvin(temperature[i_exp]);
          
    if(RH < ERH) {
      log_concentration_factor_predicted[i_exp] =
        predict_log_concentration_factor(ERH,
                                         log_shape_parameter_concavity,
                                         log_shape_parameter_steepness,
                                         initial_mass_fraction_solute);

      log_eq_mass_predicted[i_exp] =
        predict_log_equilibrium_mass(ERH,
                                     log_shape_parameter_concavity,
                                     log_shape_parameter_steepness,
                                     initial_mass[i_exp],
                                     initial_mass_fraction_solute);
      
      decay_rate[i_exp] = exp(log_arrhenius_rate(log_A_dry,
                                                 E_a_dry,
                                                 T));
    } else {
      real log_eq_conc;
      
      log_eq_mass_predicted[i_exp] =
        predict_log_equilibrium_mass(RH,
                                     log_shape_parameter_concavity,
                                     log_shape_parameter_steepness,
                                     initial_mass[i_exp],
                                     initial_mass_fraction_solute);

      log_concentration_factor_predicted[i_exp] =
        predict_log_concentration_factor(rel_humidity[i_exp],
                                         log_shape_parameter_concavity,
                                         log_shape_parameter_steepness,
                                         initial_mass_fraction_solute);
      
      if(empirical_salt){ 
        log_eq_conc = log_concentration_factor[i_exp];
      } else {
        log_eq_conc = log_concentration_factor_predicted[i_exp];
      }

      decay_rate[i_exp] = exp(log_eq_conc +
                              log_arrhenius_rate(log_A_wet,
                                                 E_a_wet,
                                                 T));
    }

  }

