real concentration_factor_linear_evaporation(real beta,
                                             real time,
                                             real initial_water_mass){

  real current_water_mass = initial_water_mass - beta * time;
  return initial_water_mass / current_water_mass;
}
