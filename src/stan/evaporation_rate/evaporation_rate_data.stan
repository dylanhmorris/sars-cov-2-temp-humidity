  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 0> n_evap_datapoints;
  int<lower = 0, upper = n_evap_datapoints> n_used_evap_datapoints;
  int<lower = 1> n_evap_classes;
  int<lower = 1, upper = n_evap_classes>
    evap_class_id[n_evap_datapoints];
  int<lower = 1, upper = n_evap_classes> experiment_evap_class[n_experiments];
  vector<lower = 0>[n_evap_datapoints] evap_time;
  vector<lower = 0>[n_evap_datapoints] measured_mass;
  vector<lower = 0>[n_evap_classes] equilibrium_mass;
  vector<lower = 0>[n_evap_classes] initial_mass;
  ////////////////////////////
  // prior hyperparameters
  ////////////////////////////
  real mean_drying_time;
  real<lower = 0> sd_drying_time;

  real mode_sd_evap;
  real<lower = 0> sd_sd_evap;

  //flags
  int <lower = 0, upper = 1> mass_prior_check;
