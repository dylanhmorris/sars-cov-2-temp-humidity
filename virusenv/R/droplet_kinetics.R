##################################################
## filename: droplet_kinetics.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
##
## description: functions for analyzing the
## kinetics of droplet drying, according
## to mechanistic chemical models
##################################################

##################################################
## functions
##################################################

#' pitzer_phi
#'
#' Solve the Pitzer equation for the
#' practical osmotic coefficient Phi at
#' a given temperature and molality for
#' a 1-1 solute
#' 
#' @param molality molality of solute (moles / kg solvent)
#' @param T_kelvin ambient temperature in Kelvin
#' @param B_0 Pitzer B_0 parameter
#' @param B_1 Pitzer B_1 parameter
#' @param C Pitzer C_parameter
#' @param pressure_atm ambient pressure in atmospheres, default 1
#' @return practical osmotic coefficient Phi
#' @export
pitzer_phi <- function(molality,
                       T_kelvin,
                       B_0,
                       B_1,
                       C,
                       pressure_atm = 1){
    alpha_1 <- 2
    I <- molality
    A_phi <- A_phi(T_kelvin, pressure_atm)
    sqrtI <- sqrt(I)
    f_phi <- A_phi * sqrtI / (1 + 1.2 * sqrtI)
    B_phi <- B_0 + B_1 * exp(-alpha_1 * sqrtI)
    Z <- 2 * molality
    C_val <- C / (2)
    return( 1 -
            f_phi +
            molality * (B_phi + Z * C_val) )
}


#' pitzer_phi_nacl
#'
#' Solve the Pitzer equation for the
#' practical osmotic coefficient Phi of
#' NaCl at a given temperature and molality
#' @param molality molality of NaCl (moles / kg solvent)
#' @param T_kelvin ambient temperature in Kelvin
#' @param pressure_atm ambient pressure in atmospheres, default 1
#' @return practical osmotic coefficient Phi
#' @export
pitzer_phi_nacl <- function(molality,
                            T_kelvin,
                            pressure_atm = 1){
    ## parameters per
    ## https://pubs.acs.org/doi/pdf/10.1021/acs.jced.6b00236

    return(pitzer_phi(molality,
                      T_kelvin,
                      B_0_NaCl,
                      B_1_NaCl,
                      C_NaCl,
                      pressure_atm))
}

#' density_pure_water
#'
#' Get the predicted density of pure water
#' at a given temperature and pressure,
#' using equations from PHREEQC
#' @param T_kelvin ambient temperature in Kelvin
#' @param pressure_atm ambient pressure in atmospheres, default 1
#' @export
density_pure_water <- function(T_kelvin,
                               pressure_atm = 1){
    ## Calculate approx density of pure water
    ## at the given T and P

    ## follows public domain C++ code
    ## from PHREEQC
    
    ## Wagner and Pruss, 2002, JPCRD 31, 387, eqn. 2.6, along the saturation pressure line +
    ## interpolation 0 - 300 oC, 0.006 - 1000 atm...
    T <- T_kelvin
    tc <- T_kelvin - 273.15
    stopifnot( tc < 350 )

    pa <- pressure_atm
    Tc <- 647.096;
    th <- 1 - T / Tc
    b1 <- 1.99274064
    b2 <- 1.09965342
    b3 <- -0.510839303
    b4 <- -1.75493479
    b5 <- -45.5170352
    b6 <- -6.7469445e5

    rho_0_sat <- (322.0 * (1.0 +
                           b1 * th^(1/3) +
                           b2 * th^(2/3) +
                           b3 * th^(5/3) +
                           b4 * th^(16/3) +
                           b5 * th^(43/3) +
                           b6 * th^(110/3)))
        
    ## pressure
    p0 <- (5.1880000E-2 +
           tc * (-4.1885519e-4 +
                 tc * ( 6.6780748e-6 + tc * (-3.6648699e-8 +
                                             tc * 8.3501912e-11))))
    
    p1 <- (-6.0251348E-06 + tc * ( 3.6696407E-07 + tc * (-9.2056269E-09 + tc * ( 6.7024182E-11 + tc * -1.5947241E-13))))
    p2 <- (-2.2983596E-09 + tc * (-4.0133819E-10 + tc * ( 1.2619821E-11 + tc * (-9.8952363E-14 + tc *  2.3363281E-16))))
    p3 <- (7.0517647E-11 + tc * ( 6.8566831E-12 + tc * (-2.2829750E-13 + tc * ( 1.8113313E-15 + tc * -4.2475324E-18))))
    
    p_sat <- exp(11.6702 - 3816.44 / (T - 46.13))

    pa <- pa - (p_sat - 1e-6)

    rho_0 <- rho_0_sat + pa * (p0 + pa * (p1 + pa * (p2 + sqrt(pa) * p3)));

    stopifnot( rho_0 > 0.01 )

    return (rho_0 / 1e3)
}


#' A_phi
#'
#' Get the predicted A_phi parameter for the Pitzer equation
#' at a given temperature and pressure,
#' using equations from PHREEQC
#' @param T_kelvin ambient temperature in Kelvin
#' @param pressure_atm ambient pressure in atmospheres, default 1
#' @export
A_phi <- function(T_kelvin,
                  pressure_atm = 1){
    ## follows public domain code from
    ## PHREEQC

    ## declare empirical constants
    n_avogadro <- 6.02252e23
    u1 <- 3.4279e2
    u2 <- -5.0866e-3
    u3 <- 9.469e-7
    u4 <- -2.0525
    u5 <- 3.1159e3
    u6 <- -1.8289e2
    u7 <- -8.0325e3
    u8 <- 4.2142e6
    u9 <- 2.1417

    ## shorten variable names
    T <- T_kelvin
    pa <- pressure_atm
    
    ## relative dielectric constant at 1000 bar
    d1000 <- u1 * exp(T_kelvin * (u2 + T * u3))
    
    c <- u4 + u5 / (u6 + T)
    b <- u7 + u8 / T + u9 * T
    pb <- pa * 1.01325 ## pa in bar
    
    ## relative dielectric constant at given pressure
    eps_r <- d1000 + c * log((b + pb) / (b + 1e3)) 

    ## constant term incorporating fundamental charge
    ## etc
    e2_DkT <- 1.671008e-3 / (eps_r * T)

    ## density of pure water
    rho_0 <- density_pure_water(T, pa)

    ##  Debye length parameter, 1/cm(mol/kg)^-0.5
    DH_B <- sqrt(8 * pi * n_avogadro * e2_DkT * rho_0 / 1e3)
    A_phi <- DH_B * e2_DkT / 6.0

    return(A_phi)
}

#' log_pred_molality_pitzer
#'
#' Get approximate log equililbrium salt molality (kg/mol)
#' as a function of relative humidity and salt constant
#' product v (num ions) * phi (practical osmotic coefficient)
#'#'
#' @param RH relative humidity
#' @param phi_salt practical osmotic coefficient
#' @param v_salt number of ions, default 2 (NaCl)
#' @param molar_mass_solvent molar mass (kg/mol) of the solvent, default 0.0180153 (H20)
#' @return log of the predicted equilibrium salt concencration
#' @export
log_pred_molality_pitzer <- function(RH,
                              phi_salt,
                              v_salt = 2,
                              molar_mass_solvent = 0.0180153){
    divisor <- v_salt * phi_salt * molar_mass_solvent
    return( log(log(RH) / (-divisor)) )
}

#' log_eq_molality_pitzer
#'
#' Solve for the log equililbrium salt molality (kg/mol)
#' as a function of relative humidity by iterating the
#' \code{\link{log_pred_molality}} function with a practical
#' osmotic coefficient phi calculated per a simplified
#' Pitzer equation to find a fixed point
#' @param RH relative humidity
#' @param v_salt number of ions, default 2 (NaCl)
#' @param molality_init_guess initial guess for the answer, default 1
#' @param molar_mass_solvent molar mass (kg/mol) of the solvent, default 0.0180153 (H20)
#' @param pressure_atm ambient pressure in atmosphers, default 1
#' @param B_0 Pitzer B_0 parameter, default 0.07722 (NaCl)
#' @param B_1 Pitzer B_1 parameter, default 0.25183 (NaCl)
#' @param C Pitzer C_parameter, default 0.00106 (NaCl)
#' @param max_iter maximum number of iterations to use, default 1000
#' @param tolerance tolerance for the log estimate, default 10^(-2)
#' @return solution for the log equilibrium salt concencration
#' @export
log_eq_molality_pitzer <- function(RH,
                                   T_kelvin,
                                   v_salt = 2,
                                   molality_init_guess = 1,
                                   molar_mass_solvent = 0.0180153,
                                   pressure_atm = 1,
                                   B_0 = B_0_NaCl,
                                   B_1 = B_1_NaCl,
                                   C = C_NaCl,
                                   max_iter = 1000,
                                   tolerance = 1e-2){
    
    guess_old <- molality_init_guess
    iter <- 1
    guess_diff <- 100 * tolerance
    
    while(
    (guess_diff > tolerance) & (iter < max_iter) ) {
        guess_phi <- pitzer_phi(guess_old,
                                T_kelvin,
                                B_0 = B_0,
                                B_1 = B_1,
                                C = C,
                                pressure_atm = pressure_atm)
        
        guess_new <- log_pred_molality_pitzer(RH,
                                              guess_phi,
                                              v_salt,
                                              molar_mass_solvent)
        guess_diff <- abs(guess_old - guess_new)
        iter <- iter + 1

        guess_old <- guess_new
    }

    stopifnot(iter < max_iter)

    return (guess_old)
}


#' molality
#' 
#' Get the molality of a solute (default NaCl).
#' Makes the assumption that the volume occupied
#' by the solute can be neglected
#' 
#' @param solute_grams mass of solute in grams
#' @param volume_liters volume of solvent in liters
#' @param solvent_density_kg_L density of the solvent in kg per liter
#' @param solute_molar_mass_grams molar mass of the solute in g/mol
#' @return molality of the solution
#' @export
molality <- function(solute_grams,
                     volume_liters = 1,
                     solvent_density_kg_L = 1, # H20 STP
                     solute_molar_mass_grams = 58.44 # NaCl, g/mol
                     ){
    moles_per_vol <- solute_grams / solute_molar_mass_grams
    solvent_mass_kg <- (volume_liters * solvent_density_kg_L -
                        solute_grams / 1000)
    return(moles_per_vol / solvent_mass_kg)
}

#' mass_frac
#' 
#' Get the mass fraction of a solute (default NaCl)
#' Makes the assumption that the volume occupied
#' by the solute can be neglected
#' 
#' @param solute_grams mass of solute in grams
#' @param volume_liters volume of solvent in liters
#' @param solvent_density_kg_L density of the solvent in kg per liter
#' @return mass fraction of the solute
#' @export
mass_frac <- function(solute_grams,
                      volume_liters = 1,
                      solvent_density_kg_L = 1 # H20 STP
                      ){
    solute_mass_kg <- (solute_grams / 1000)
    total_mass_kg <- (volume_liters * solvent_density_kg_L +
                      solute_mass_kg)
    return(solute_mass_kg / total_mass_kg)
}

#' mass_frac_to_molality
#' 
#' Get the molality of a solute (default NaCl in water)
#' given the mass fraction
#' 
#' @param mass_fract mass fraction of solute
#' @param solute_molar_mass_grams molar mass of the solute in
#' grams (default 58.44 (NaCl))
#' @return molality of the solute
#' @export
mass_frac_to_molality <- function(mass_frac,
                                  solute_molar_mass_grams = 58.44 # NaCl
                                  ){
    solute_mass_kg <- mass_frac / (1 - mass_frac)
    solute_mass_g <- solute_mass_kg * 1000
    solute_moles <- solute_mass_g / solute_molar_mass_grams

    return(solute_moles)
}


#' molality_to_mass_frac
#' 
#' Get the mass fraction of a solute (default NaCl in water)
#' given the molality
#' 
#' @param molality molality of solute
#' @param solute_molar_mass_grams molar mass of the solute in
#' grams (default 58.44 (NaCl))
#' @return mass fraction of the solute
#' @export
molality_to_mass_frac <- function(molality,
                                  solute_molar_mass_grams = 58.44 # NaCl
                                  ){
    solute_mass_kg <- molality * solute_molar_mass_grams / 1000

    return(solute_mass_kg / (1 + solute_mass_kg))
}


#' cohen_RH_of_molality
#' 
#' Get the relative humidity
#' corresponding to a given molality of
#' NaCl in pure water at env conditions
#' from an empirical equation per Cohen et al
#' https://doi.org/10.1021/j100301a029
#' @param molality molality of NaCl
#' @return corresponding relative humidity
#' @export
cohen_RH_of_molality <- function(molality){
    m <- molality
    return(1.0084 -
           4.939e-2 * m +
           8.888e-3 * m^2 -
           2.157e-3 * m^3 +
           1.617e-4 * m^4 -
           1.990e-6 * m^5 -
           1.142e-7 * m^6)
}

#' cohen_molality_of_RH
#' 
#' Get the molality of
#' NaCl in pure water at env conditions
#' corresponding to the given RH by numerically
#' inverting the empirical equation per Cohen et al
#' https://doi.org/10.1021/j100301a029
#' @param RH ambient relative humidity
#' @return corresponding molality of NaCl
#' @export
cohen_molality_of_RH <- function(RH,
                           tolerance = 1e-2){
    result <- uniroot(
        function(x){cohen_RH_of_molality(x) - RH},
        c(0, 2000))

    stopifnot(result$estim.prec < tolerance)

    return (result$root)
}

#' tang_aw
#' 
#' Get the water activity of
#' NaCl in pure water at env conditions
#' corresponding to the given mass fraction
#' using the empirical equation of Tang 1996
#' https://doi.org/10.1029/96JD03003
#' @param mass_fract mass fraction of NaCl
#' @return corresponding water activity
#' @export
tang_aw <- function(mass_frac){
    ## Tang uses percentages, not decimals
    mass_frac <- mass_frac * 100
    vars_vec <- c(1,
                  mass_frac,
                  mass_frac^2,
                  mass_frac^3,
                  mass_frac^4)
    return( sum(tang_NaCl_vec * vars_vec) )
}


#' tang_mass_frac
#' 
#' Get the mass fraction of NaCl
#' NaCl in pure water at the given RH
#' by inverting the empirical equation of Tang 1996
#' https://doi.org/10.1029/96JD03003
#' analytically using the quartic formula
#' 
#' @param RH ambient relative humidity, as a decimal
#' @return mass fraction of NaCl
#' @export
tang_mass_frac <- function(RH){
    if((RH < tang_bounds_NaCl[1]) |
       (RH > tang_bounds_NaCl[2])){
        stop("RH out of bounds for Tang equation")
    }
    
    vec_to_solve <- tang_NaCl_vec - c(RH, rep(0, 4))
    
    return( solve_tang_quartic(vec_to_solve) / 100)
}


#' mass_frac_to_molar_frac
#'
#' convert a mass fraction to a corresponding
#' molar fraction (default NaCl in water)
#'  
#' @param mass_frac mass fraction of the component of interests
#' @param molar_mass molar mass of the component of interest,
#' default 0.05844 mol/kg (NaCl)
#' @param molar_mass_others average molar mass of other components,
#' default 0.0180153 mol/kg (pure H20)
#' @export
mass_frac_to_molar_frac <- function(mass_frac, # H20 STP
                                    molar_mass = 0.05844,    #NaCl
                                    molar_mass_others = 0.0180153) # pure H20
{
    moles <- mass_frac  * molar_mass
    moles_other <- molar_mass_others * (1 - mass_frac)

    return( moles / (moles + moles_other) )
}

#' fraction_change_to_log_concentration_factor
#'
#' convert a change in mass fraction or molar
#' fraction to a log concentration factor (log fold-concentration
#' relative to initial)
#'  
#' @param fraction_final final fraction of the component of interest
#' @param fraction_initial initial fraction of the component of interest
#' @return natural log of how many times concentrated that quantity has been
#' @export
fraction_change_to_log_concentration_factor <- function(fraction_final,
                                                        fraction_initial){
    init_term <- log(1 - fraction_initial) - log(fraction_initial)
    final_term <- log(fraction_final) - log(1 - fraction_final)
    return( init_term + final_term )
}

#' mass_change_to_log_concentration_factor
#'
#' convert a change in mass to a log
#' concentration factor (log fold-concentration
#' relative to initial)
#'  
#' @param mass_final final mass of the solution
#' @param mass_initial initial mass of the solution
#' @param fraction initial initial mass fraction of the
#' conserved solution components (solutes, usually)
#' @return natural log of how many times concentrated that quantity has been
#' @export
mass_change_to_log_concentration_factor <- function(mass_final,
                                                    mass_initial,
                                                    fraction_initial){
    conserved_mass <- fraction_initial * mass_initial

    return( log(mass_initial - conserved_mass) -
            log(mass_final - conserved_mass))
}


#' solute_to_solvent_mole_ratio
#'
#' predict the ratio of solute mole
#' fraction to solvent mole fraction
#' at evaporative equilibrium for a
#' given relative humidity, using a
#' two-parameter semi-mechanistic
#' functional form
#'  
#' @param relative_humidity ambient relative humidity
#' @param log_shape_parameter_concavity parameter governing
#' concavity of the mole-fraction / water activity relationship (negative
#' = concave down, positive = concave up)
#' log_shape_parameter_steepness parameter governing
#' steepness of the mole-fraction / water activity relationship near a_w/rh = 1(larger values = steeper)
#' @return predicted ratio of solute mole fraction to
#' solvent mole fraction at equilibrium
#' @export
solute_to_solvent_mole_ratio <- function(relative_humidity,
                                         log_shape_parameter_concavity,
                                         log_shape_parameter_steepness){
    a_conc <- exp(-log_shape_parameter_concavity)
    inv_conc <- 1 / a_conc
    a_steep <- exp(-log_shape_parameter_steepness)
    
    return( (-log(relative_humidity) / a_steep)^inv_conc )
}


#' predict_log_concentration_factor
#'
#' predict (the natural log of)
#' how many-fold concentrated a solution
#' will be at equilibrium given
#' the initial mass fraction,
#' two shape parameters governing
#' non-ideality
#' and the initial solute concentration
#'  
#' @param relative_humidity ambient relative humidity
#' @param log_shape_parameter_concavity parameter governing
#' concavity of the mole-fraction / water activity relationship (negative
#' = concave down, positive = concave up)
#' log_shape_parameter_steepness parameter governing
#' steepness of the mole-fraction / water activity relationship (larger values = steeper)
#' @param initial_mass_fraction_solute mass fraction of the
#' solute(s) in the solution at t = 0.
#' @return log concentration factor for the solution
#' @export
predict_log_concentration_factor <- function(relative_humidity,
                                             log_shape_parameter_concavity,
                                             log_shape_parameter_steepness,
                                             initial_mass_fraction_solute){
  
    log_final_ratio <- (
        log(solute_to_solvent_mole_ratio(relative_humidity,
                                         log_shape_parameter_concavity,
                                         log_shape_parameter_steepness)))

    log_initial_ratio <- (log(initial_mass_fraction_solute) -
                          log(1 - initial_mass_fraction_solute))
    
    return(log_final_ratio - log_initial_ratio)
}



#' buck_eqn
#' 
#' returns a value for saturation vapor
#' pressure in Pa using Buck's empiricial
#' equation
#' 
#' @param T_kelvin Ambient temperature in Kelvin
#' @return corresponding saturation vapor pressure
#' @export
buck_eqn <- function(T_kelvin){
    return( 6.1121e-2 * exp((18.678 - T_kelvin / 234.5) *
                            (T_kelvin / (257.14 + T_kelvin))))
}

#' to_AH
#' 
#' convert a value of RH to its AH given
#' the ambient temperature T_kelvin 
#'
#' @param RH relative humidity
#' @param T_kelvin Ambient temperature in Kelvin
#' @return corresponding absolute humidity in g/m^3
#' @export
to_AH <- function(RH, T_kelvin){
    C <- 2.16679 # gK/J
    P_sat <- buck_eqn(T_kelvin)
    P_w <- P_sat * RH
    AH <- C * P_w / T_kelvin
}
