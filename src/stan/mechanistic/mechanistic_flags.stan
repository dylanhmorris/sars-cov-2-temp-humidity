  // whether to use measured or modeled
  // salt concentrations for RH > ERH
  int<lower = 0, upper = 1> empirical_salt;

  // which function form to use for model
  // concentration
  int<lower = 1> functional_form;
