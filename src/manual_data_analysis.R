library(rstan)
library(shinystan)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(virusenv)
library(Rcpp)
library(tidyr)
library(readr)
library(viridis)


qpcr <- "../dat/cleaned/aerosol_metadata.csv"

current <- '../out/mcmc_chains/mechanistic-evaporation-prior-check.Rds'



current <- "../out/mcmc_chains/phi6_chains.Rds"
current <- "../out/mcmc_chains/aerosol_prior_check.Rds"
current <- "../out/mcmc_chains/mechanistic_plastic_chains.Rds"
current <- "../../out/mcmc_chains/mechanistic_evaporation_chains.Rds"
current <- "../out/mcmc_chains/inferred_aerosol_titers_decay_chains.Rds"
current <- "../out/mcmc_chains/aerosol_humidity_bin_chains.Rds"
pla_output_path <- "../../covid19/out/mcmc_chains/th_plastic_model_salt_chains.Rds"
pla_hl_path <- '../../out/mcmc_chains/plastic_chains.Rds'

lit_output_path <- "../out/mcmc_chains/literature_basic_chains.Rds"
aer_inf_output_path <- "../out/mcmc_chains/inferred_aerosol_titers_decay_chains.Rds"
aer_inf_output_path <- "../out/mcmc_chains/inferred_aerosol_titers_decay_chains.Rds"

p_check_path <- "../out/mcmc_chains/inferred_aerosol_titers_decay_chains.Rds"

aer_dat <- read_csv("../dat/cleaned/th_aerosol_data.csv")
a_chains <- readRDS(aer_inf_output_path)

hlh <- function(reference_half_life,
                temp_pred,
                temp_ref,
                hum_pred,
                hum_ref,
                ERH,
                E_a_dry,
                E_a_wet,
                log_concavity_parameter = 0){

    super <- hum_ref > ERH

    inv_a_conc <- 1 / exp(log_concavity_parameter)

    conc_ratio <- ifelse(
        super,
        (log(hum_ref) / log(hum_pred))^(inv_a_conc),
        1)

    E_a <- ifelse(
        super,
        E_a_wet,
        E_a_dry)
    
    temp_inv_diff <- (1 / temp_pred) - (1 / temp_ref)

    temp_ratio <- exp(E_a * temp_inv_diff / R)
    
    pred_hl <- reference_half_life * conc_ratio * temp_ratio

    return (pred_hl)
}



c %>% ggplot(aes(
          x = time,
          y = 10^(log10_fraction_remaining),
          group = titer_id)) +
    geom_violin() +
    facet_grid(vars(target_humidity,
                    replicate),
               vars(target_temperature)) +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(ylim = c(1/500, 2))


aer_titers <- aer_dat %>%
    distinct(titer_id, .keep_all = TRUE) %>%
    mutate(AH = to_AH(humidity_measured,
                      temperature_measured + 273.15)) %>%
    arrange(time) %>%
    group_by(run_number, .draw) %>%
    mutate(
        previous_AH = lag(AH,
                          default = NA),
        log_humidity_drop = log(previous_humidity) - log(humidity_measured)
    )


a_draws <- a_chains %>% spread_draws(true_titer[titer_id]) %>%
    inner_join(aer_dat %>% distinct(titer_id, .keep_all = TRUE),
               by = "titer_id")

a_draws %>%
    ggplot(aes(
        x = time,
        y = 10^(true_titer + 1),
        fill = humidity_measured,
        group = run_number)) +
    stat_eye(width = 2,
             scale = 0.5) +
    scale_fill_viridis() +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(ylim = c(1, 5e3)) + 
    facet_grid(vars(target_temperature),
               vars(target_humidity, replicate))

chains <- readRDS(pla_output_path) %>%
    spread_draws(log_init_mass_frac[experiment_id])

ERH <- 0.45

a = read_csv("../../covid19/dat/cleaned/th_plastic_data.csv") %>%
    distinct(trial_unique_id, .keep_all = TRUE) %>%
    crossing(chains) %>%
    mutate(
        log_concentration_factor = log(initial_mass) - log(equilibrium_mass),
        molality_initial = log_init_mass_frac %>% exp() %>% mass_frac_to_molality(),
        molality_final = exp(log_concentration_factor + log_init_mass_frac) %>% mass_frac_to_molality()) %>%
    select(trial_unique_id,
           humidity,
           temperature,
           log_concentration_factor,
           molality_initial,
           molality_final,
           .draw)

b = tibble(
    rh = rep(seq(0.4, 0.999, length.out = 100), 3)) %>%
    crossing(tibble(temperature = c(10, 22, 27))) %>%
    rowwise() %>%
    mutate(molality_mech = (
        exp(ifelse(rh > ERH,
                   log(20.4245628903 * sqrt(-log(rh))),
                   log_eq_molality(ERH,
                                   temperature + 273.1)))) %>%
            molality_to_mass_frac() %>%
            molarity_from_mass_frac(),
        
        molality_marr = (
            ifelse(rh > ERH,
                   cohen_molality_of_RH(rh),
                   cohen_molality_of_RH(ERH)) %>%
            molality_to_mass_frac() %>%
            molarity_from_mass_frac()),

        molality_tang = ifelse(rh > ERH,
                               tang_mass_frac(rh),
                               tang_mass_frac(ERH)) %>%
            molarity_from_mass_frac())


a %>% filter(.draw %in% sample(1:4000, 500)) %>%
    ggplot(
        aes(x = humidity / 100,
            y = molality_final,
            color = temperature,
            fill = temperature)) +
    geom_point(alpha = 0.05, position = "jitter") +
    geom_point(aes(x = 1,
                   y = molality_initial),
               alpha = 0.05, position = "jitter") +
    scale_fill_viridis() + 
    scale_color_viridis() +
    geom_line(data = b,
              mapping = aes(
                  x = rh,
                  y = molality_mech),
              color = "black") +
    geom_line(data = b,
              mapping = aes(
                  x = rh,
                  y = molality_marr),
              color = "black",
              linetype = "dashed") +
    geom_line(data = b,
              mapping = aes(
                  x = rh,
                  y = molality_tang),
              color = "black",
              linetype = "dotted") +    
    geom_pointinterval() +
    scale_color_viridis() +
    scale_y_continuous(trans = "log2") +
    coord_cartesian(xlim = c(0, 1))



a = readRDS(current) %>%
    spread_draws(half_life[th_bin_id]) %>%
    add_titer_metadata(dat, "th_bin_id") %>%
    mutate(prior_density = dnorm(log(half_life), log(3), 1.25))

a %>% ggplot(aes(x = half_life, y = humidity_bin, pointfill = target_temperature)) + facet_wrap(vars(target_temperature)) + stat_halfeye(shape = 21, quantiles = 100) + scale_fill_viridis() + geom_line(aes(y = prior_density)) + coord_cartesian(xlim = c(0, 30))
