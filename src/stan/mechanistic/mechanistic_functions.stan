/* filename: mechanistic_model.stan
 * author: Dylan Morris <dhmorris@princeton.edu>
 * description: defines mechanistic model
 * of droplet evaporation and virus inactivation
 * used for fitting multiple stan models
 */


/////////////////////////////////
// physical constants
/////////////////////////////////
  
   // gas constant in J / Kelvin mol
real R(){
  return 8.31446261815324;
}

// 0 C in Kelvin
real zero_C_in_Kelvin(){
  return 273.15;
}

real to_Kelvin(real temp_C){
  return temp_C + zero_C_in_Kelvin();
}
  
real to_Celsius(real temp_K){
  return temp_K - zero_C_in_Kelvin();
}

/////////////////////////////////
// convenience math
/////////////////////////////////

// difference between 1/A and 1/B in
// kelvin where A and B are in Celsius
real temp_inv_diff(real temp_one,
                   real temp_two){
  real one = temp_one + zero_C_in_Kelvin();
  real two = temp_two + zero_C_in_Kelvin();
  return (1 / one - 1 / two);
}

  
real tang_mass_frac (real RH){
  real e_t = 1 - RH;
  real d = (-6.366e-3);
  real c = (8.624e-5);
  real b = (-1.158e-5);
  real a = (1.518e-7);    

  real p = (8.0 * a * c - 3.0 * pow(b, 2)) / (8.0 * pow(a, 2));
  real q = (pow(b, 3) - 4.0 * a * b * c + 8.0 * pow(a, 2) * d) /
    (8.0 * pow(a, 3));
  
  real Delta0 = pow(c, 2) - 3.0 * b * d + 12 * a * e_t;
  real Delta1 = (2 * pow(c, 3) -
                 9.0 * b * c * d +
                 27.0 * pow(b, 2) * e_t +
                 27.0 * a * pow(d, 2) -
                 72.0 * a * c * e_t);

  real Q =
    pow(((Delta1 + sqrt(pow(Delta1, 2) -
                        4.0 * pow(Delta0, 3))) / 2.0),
        1.0 / 3.0);
  
  real S = 0.5 * sqrt((-2.0 * p / 3.0) +
                      (1.0 / (3.0 * a)) *
                      (Q + Delta0 / Q));
  
  real term1 = -b / (4.0 * a);
  real term3 = 0.5 * sqrt(-4.0 * pow(S, 2) - 2.0 * p - q / S);
  // one of the roots is always real and our desired solution
  // namely, this one 
  real x4 = term1 + S - term3;
  
  // tang operates in percentages; we operate in decimals
  return x4 / 100;
}


// given ratio x/(1-x) where 0 < x < 1, find x
real ratio_to_unit_fraction(real ratio){
  return ratio / (ratio + 1);
}


real solute_to_solvent_mole_ratio(real relative_humidity,
                                  real log_shape_parameter_concavity,
                                  real log_shape_parameter_steepness){

  real a_conc = exp(-log_shape_parameter_concavity);
  real inv_conc = 1 / a_conc;
  
  real a_steep = exp(-log_shape_parameter_steepness);
  
  return pow(-log(relative_humidity) / a_steep, inv_conc);
}


real predict_log_concentration_factor(real relative_humidity,
                                      real log_shape_parameter_concavity,
                                      real log_shape_parameter_steepness,
                                      real initial_mass_fraction_solute){
  
  real log_final_ratio = log(
    solute_to_solvent_mole_ratio(relative_humidity,
                                 log_shape_parameter_concavity,
                                 log_shape_parameter_steepness));
  real log_initial_ratio = log(initial_mass_fraction_solute) -
    log(1 - initial_mass_fraction_solute);
  
  return log_final_ratio - log_initial_ratio;
}


real predict_log_equilibrium_mass(real relative_humidity,
                                  real log_shape_parameter_concavity,
                                  real log_shape_parameter_steepness,
                                  real initial_mass,
                                  real initial_mass_fraction_solute){
  real solute_mass_initial = initial_mass * initial_mass_fraction_solute;
  real fraction_final = ratio_to_unit_fraction(
    solute_to_solvent_mole_ratio(relative_humidity,
                                 log_shape_parameter_concavity,
                                 log_shape_parameter_steepness));
  return log(solute_mass_initial) - log(fraction_final);
}

/////////////////////////////////
// mechanistic model 
/////////////////////////////////

real log_arrhenius_rate(real log_A,
                        real E_a,
                        real T_kelvin){
  return log_A - E_a / (R() * T_kelvin);
}

real log_arrenhius_A(real k_standard,
                     real E_a,
                     real T_standard_kelvin){
  return log(k_standard) + E_a / (R() * T_standard_kelvin);
}



