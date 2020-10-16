  for(i_exp in 1:n_experiments) {
    real RH = rel_humidity[i_exp];
    real T = to_Kelvin(temperature[i_exp]);
    
    transient_decay_rate[i_exp] =
      exp(log_arrhenius_rate(log_A_wet,
                             E_a_wet,
                             T));

  }

