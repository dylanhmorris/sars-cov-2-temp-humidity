  class_decay_rate = log10(2) ./ exp(transient_log_half_life);

  for(i_exp in 1:n_experiments){
    int i_class = transient_class_id[i_exp];
    transient_decay_rate[i_exp] = class_decay_rate[i_class];
  }

