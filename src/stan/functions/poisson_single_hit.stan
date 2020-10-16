  /* filename: poisson_single_hit.stan
   * author: Dylan H. Morris <dhmorris@princeton.edu>
   * 
   * description: header file defining function
   * to increment the log probability according to
   * a poisson single hit model
   */

real poisson_single_hit_lpmf(int well_status,
                             real virions){
  real result = 0;
  
  if (well_status == 0) {
    result = poisson_lpmf(0 | virions);
  } else if (well_status == 1) {
    result = poisson_lccdf(0 | virions);
  } else {
    reject("well_status data must be one or zero; given",
           well_status);
  }
  
  return result;
  
}
