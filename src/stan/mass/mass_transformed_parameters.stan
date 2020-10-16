real initial_mass_fraction_solute = exp(log_initial_mass_solute_per_thousand) / 1000;

vector[n_evap_classes] solute_mass_init = initial_mass_fraction_solute * initial_mass;

vector[n_evap_classes] log_concentration_factor = log(initial_mass - solute_mass_init) - log(equilibrium_mass - solute_mass_init);
