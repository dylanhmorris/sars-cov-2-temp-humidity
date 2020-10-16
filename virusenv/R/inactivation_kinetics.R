##################################################
## filename: inactivation_kinetics.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
##
## description: functions for analyzing the
## kinetics of virus inactivation, according
## to mechanistic chemical models
##################################################

##################################################
## functions
##################################################


#' log_arrhenius_rate
#'
#' Get the log reaction rate from the Arrenhius equation
#' given log asymptotic reaction rate log(A), activation energy
#' E_a and temperature in kelvin
#' @param log_A log of the asymptotic reaction rate at high temperature
#' @param E_a activation energy for the reaction in J/(Kelvin mol)
#' @param T_kelvin temperature in Kelvin
#' @return log of the predicted equilibrium salt concentration
#' @export
log_arrhenius_rate <- function(log_A,
                               E_a,
                               T_kelvin){
    return ( log_A - E_a / (R * T_kelvin) )
}


#' log_arrenhius_A
#'
#' Get the log asymptotic reaction rate log(A) rate
#' given the reaction rate at some standard temperature,
#' the activation energy E_a, and the value in Kelvin of
#' the standard temperature
#' @param k_standard reaction rate at the standard temperature
#' @param E_a activation energy for the reaction in J/(Kelvin mol)
#' @param T_kelvin value of the standard temperature in Kelvin
#' @return log of the predicted equilibrium salt concencration
#' @export
log_arrenhius_A <- function(k_standard,
                            E_a,
                            T_standard_kelvin){
    return ( log(k_standard) + E_a / (R * T_standard_kelvin) )
}



#' mechanistic_decay_rate
#' 
#' Get the predicted decay rate (reaction rate of
#' the inactivation reaction) of viable virus
#' given temperature and humidity and the
#' ERH and mass fraction of the inactivating
#' reactant (salt or something else).
#' By default, mass fractions are converted
#' into molarities according to the function
#' \link{mass_frac_to_molarity}
#' 
#' @param relative_humidity ambient relative humidity as a decimal
#' @param T_kelvin ambient temperature in Kelvin
#' @param log_concentration_term log of the concentration term modifying the reaction rate
#' @param ERH efflorescence relative humidity for the inactivating agent(s)
#' @param log_A_dry log of the asymptotic reaction rate at high temperatures but below the ERH (i.e. when the agent(s) that are leading order virus inactivators in solution effloresce)
#' @param E_a_dry activation energy of the reaction below the ERH (i.e. when the agent(s) that are leading order virus inactivators in solution effloresce)
#' @param log_A_wet log of the asymptotic reaction rate at high temperatures above the ERH (i.e. in solution) with asymptotically high concentration of inactivating agents(s)
#' @param E_a_wet activation energy of the reaction above the ERH (i.e. in solution)
#' @return predicted exponetial decay rate of virus (reaction rate k of the inactivation reaction)
#' @export
mechanistic_decay_rate <- function(relative_humidity,
                                   T_kelvin,
                                   log_concentration_term,
                                   ERH,
                                   log_A_dry,
                                   E_a_dry,
                                   log_A_wet,
                                   E_a_wet){
    
    return(
        ifelse(
            relative_humidity < ERH,
            exp(log_arrhenius_rate(log_A_dry,
                                   E_a_dry,
                                   T_kelvin)),
            exp(log_concentration_term +
                log_arrhenius_rate(log_A_wet,
                                   E_a_wet,
                                   T_kelvin)))
    )
}

#' decay_ratio_pred
#'
#' get the predicted ratio of inactivation
#' rates for two temperatures given the
#' estimated activation energy.
#' Assumes constant humidity.
#' 
#' @param temp_pred temperature for the rate to predict, Kelvin
#' @param temp_ref temperature for the reference rate, Kelvin
#' @param E_a activation energy of the reaction
#' @return predicted ratio of the rates (predicted / reference)
#' @export
decay_ratio_pred <- function(temp_pred,
                             temp_ref,
                             E_a){
    return( exp((E_a / R) * (1 / temp_ref - 1 / temp_pred)) )
}

#' hl_ratio_pred
#'
#' get the predicted ratio of half-lives
#' for two temperatures given the
#' estimated activation energy.
#' Assumes constant humidity.
#' 
#' @param temp_pred temperature for the rate to predict, Kelvin
#' @param temp_ref temperature for the reference rate, Kelvin
#' @param E_a activation energy of the reaction
#' @return predicted ratio of the half-lifes (predicted / reference)
#' @export
hl_ratio_pred <- function(temp_pred,
                          temp_ref,
                          E_a){
    decay_ratio <- decay_ratio_pred(temp_pred,
                                    temp_ref,
                                    E_a)
    return( 1 / decay_ratio )
}


#' decay_rate_pred
#'
#' get the predicted inactivation
#' rate as a function of temperature
#' based on a known reference rate and
#' an estimated activation energy.
#' Assumes constant humidity.
#'
#' @param reference_rate known rate for a given temperature
#' @param temp_pred temperature for the rate to predict, Kelvin
#' @param temp_ref temperature for the reference rate, Kelvin
#' @param E_a activation energy of the reaction
#' @return predicted rate at the new temperature 
#' @export
decay_rate_pred <- function(reference_rate,
                            temp_pred,
                            temp_ref,
                            E_a){
    ratio <- decay_ratio_pred(temp_pred, temp_ref, E_a)
    return( reference_rate * ratio )
}

#' half_life_hum
#'
#' get the predicted half-life
#' as a function of humidity
#' and temperature based on
#' a known reference half-life and
#' activation energy, with the
#' assumption of an RH^(1/alpha_c)
#' functional form
#'
#' @param reference_half_life known half-life for a
#' given temperature/humidity
#' @param temp_pred temperature for the rate to predict, Kelvin
#' @param temp_ref temperature for the reference rate, Kelvin
#' @param hum_pred fractional humidity for the rate to predict
#' @param hum_ref fractional humidity for the reference rate
#' @param ERH fractional effloresecence relative humidity of the
#' solution containing the virions
#' @param E_a_dry activation energy of the reaction below the ERH
#' @param E_a_wet activation energy of the reaction above the ERH
#' @param log_concavity_parameter log of a parameter governing the
#' concavity of the relationship between water activity and RH at
#' RH close to 1. Negative: concave down. Positive: concave up.
#' 
#' @return predicted half-life at the new temperature/humidity
#' @export
half_life_hum <- function(reference_half_life,
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


#' half_life_pred
#'
#' get the predicted half-life
#' as a function of temperature
#' based on a known reference half-life and
#' an estimated activation energy.
#' Assumes constant humidity.
#'
#' @param reference_half_life known half-life for a given temperature
#' @param temp_pred temperature for the half-life to predict, Kelvin
#' @param temp_ref temperature for the reference half-life, Kelvin
#' @param E_a activation energy of the reaction
#' @return predicted half-life at the new temperature
#' @export
half_life_pred <- function(reference_half_life,
                           temp_pred,
                           temp_ref,
                           E_a){
    ref_decay_rate <- to_decay_rate(reference_half_life)
    ratio <- decay_ratio_pred(temp_pred, temp_ref, E_a)
    new_rate <- ref_decay_rate * ratio
    return( to_fractional_life(new_rate,
                               fraction = 1 / 2,
                               decay_rate_log_base = 10) )
}


#' to_decay_rate
#'
#' convert a an nth-life
#' (i.e. how long until the fraction in
#' question is all that remains)
#' into a decay rate in logs per hour
#' (default log base 10)
#'
#' @param fractional_life fractional life to convert
#' in a given time unit
#' @param fraction used (default 1 / 2: half-life)
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return decay_ rate in logs of the given
#' base per time unit
#' @export
to_decay_rate <- function(fractional_life,
                          fraction = 1 /2,
                          decay_rate_log_base = 10){

    log_frac <- -log(fraction) / log(decay_rate_log_base)
    return(log_frac / fractional_life)
}


#' to_fractional_life
#'
#' convert a decay rate to an nth-life
#' (i.e. how long until the fraction in
#' question is all that remains)
#'
#' @param decay_rate exponential decay rate
#' @param fraction fraction to calculate for
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return fractional life for the given fraction
#' @export
to_fractional_life <- function(decay_rate,
                               fraction,
                               decay_rate_log_base = 10){

    log_frac <- -log(fraction) / log(decay_rate_log_base)
    return(log_frac / decay_rate)
}


#' to_half_life
#'
#' convert a decay rate to a half_life
#'
#' @param decay_rate exponential decay rate
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return fractional life for the given fraction
#' @export
to_half_life <- function(decay_rate,
                         decay_rate_log_base = 10){
    fraction <- 1 / 2
    return (to_fractional_life(decay_rate,
                               fraction,
                               decay_rate_log_base))
}



#' predict_titers_exponential_decay
#'
#' predict log10 titer as a function of time under
#' the assumption of simple first-order kinetics
#' (exponential decay)
#'
#' @param time time since depsotion at which the titer is observed
#' @param intercept log10 titer at time = 0
#' @param decay_rate exponential decay rate of viable virus
#' @return predicted log10 titer
#' @export
predict_titers_exponential_decay <- function(time,
                                             intercept,
                                             decay_rate){
    return( intercept - decay_rate * time )
}


#' predict_titers_implicit_evaporation
#'
#' predict log10 titer as a function of time under
#' the assumption of two phases of
#' simple first-order kinetics
#' (exponential decay): a transient phase prior
#' to evaporative equilibrium and an equilibrium
#' phase thereafter
#'
#' @param time time since depsotion at which the titer is observed
#' @param intercept log10 titer at time = 0
#' @param transient_decay_rate exponential decay rate of viable virus
#' prior to evaporative equilibrium
#' @param eq_decay_rate exponential decay rate of viable virus at
#' evaporative equilibrium
#' @param drying_time time evaporative equilibrium is reached
#' @return predicted log10 titer
#' @export
predict_titers_implicit_evaporation <- function(time,
                                                intercept,
                                                transient_decay_rate,
                                                eq_decay_rate,
                                                drying_time){
    predicted_titers <- ifelse(
        time < drying_time,
        intercept - transient_decay_rate * time,
        (intercept -
         transient_decay_rate * drying_time -
         eq_decay_rate * (time - drying_time)))

    return (predicted_titers)
}


#' predict_titers_explicit_evaporation
#'
#' predict log10 titer as a function of time with
#' explicit modeling of evaporation during
#' the transient phase
#'
#' @param time time since depsotion at which the titer is observed
#' @param intercept log10 titer at time = 0
#' @param initial_decay_rate decay rate of viable virus at t = 0
#' @param eq_decay_rate exponential decay rate of viable virus after evaporative equilibrium is reached
#' @param normalized_evap_rate_B normalized linear rate B of evaporation, in the same units as time (so that full evaporation would occur at t = 1 / B
#' @param equilibrium_concentration_factor degree of concentration of the solution (relative to the initial) at equilibrium (sets the point at which evaporation stops)
#' @return predicted log10 titer
#' @export
predict_titers_explicit_evaporation <- function(time,
                                                intercept,
                                                initial_decay_rate,
                                                eq_decay_rate,
                                                normalized_evap_rate_B,
                                                equilibrium_concentration_factor)
{
    B <- normalized_evap_rate_B
    eq_conc <- equilibrium_concentration_factor
    t <- time
    drying_t <- (1 - (1 / eq_conc)) / B

    transient_log_term <- ifelse(
        t < drying_t,
        1 - B * t,
        1 - B * drying_t)
    
    predicted_titer <- ifelse(
        t < drying_t,
        intercept +
        (initial_decay_rate / B) *
        log10(transient_log_term),
        intercept +
        (initial_decay_rate / B) *
        log10(transient_log_term) -
        eq_decay_rate * (t - drying_t))
    
    return( predicted_titer )

}

