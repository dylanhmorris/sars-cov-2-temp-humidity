functions {
  
  real normal_ub_rng(real mu, real sigma, real ub) {
    real p_ub = normal_cdf(ub, mu, sigma);
    real u = uniform_rng(0, p_ub);
    real error = inv_Phi(u);
    real y = mu + sigma * error;
    return y;
  }
}


data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 1> n_experiments;
  int<lower = 1> n_experiment_classes;
  int<lower = 1> n_measurements;
  int<lower = 0> n_measurements_used;
  
  int<lower = 1, upper = n_experiment_classes>
    experiment_class_id[n_measurements];
  
  int<lower = 1, upper = n_experiments>
    experiment_id[n_measurements];

  vector[n_measurements] log10_fraction_measured;
  
  vector<lower = 0>[n_measurements] time;
  vector<upper = 0>[n_measurements] log10_fraction_LOD;
  
  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real log_hl_prior_mean;
  real<lower = 0> log_hl_prior_sd;
 
  real<lower = 0> mode_sd_sample;
  real<lower = 0> sd_sd_sample;
  
  int<lower = 0, upper = 1> debug;

}

transformed data {}


parameters{
  vector[n_experiments] log_half_life;
  vector<lower = 0>[n_experiment_classes] sd_sample;

}

transformed parameters {  
  vector<lower = 0>[n_experiments] decay_rate;

  decay_rate = log10(2) ./ exp(log_half_life);
}

model {

  // observation process: gaussian error in log10
  // with censoring

  for(i_meas in 1:n_measurements_used){
    int i_exp = experiment_id[i_meas];
    int i_class = experiment_class_id[i_meas];
    real log10_frac = log10_fraction_measured[i_meas];
    real mu_frac = -decay_rate[i_exp] *
      time[i_meas]; 
    real sd_frac = sd_sample[i_class];
    real LOD =  log10_fraction_LOD[i_meas];
    if (log10_frac <= LOD) {
      target += normal_lcdf((LOD - mu_frac) / sd_frac | 0, 1);
    } else {
      target +=
        std_normal_lpdf((log10_frac - mu_frac) / sd_frac);
    }
  } // close loop over observations

    
  log_half_life ~ normal(log_hl_prior_mean,
                         log_hl_prior_sd);
  
  sd_sample ~ normal(mode_sd_sample,
                     sd_sd_sample);

}

generated quantities {
  vector[n_measurements] sampled_titers_pred;

  for(i_meas in 1:n_measurements){
    int i_exp = experiment_id[i_meas];
    int i_class = experiment_class_id[i_meas];

    real ith_predicted_frac = -decay_rate[i_exp] * time[i_meas];
        
    sampled_titers_pred[i_meas] =
      normal_rng(ith_predicted_frac,
                 sd_sample[i_class]);
    
  }
}
