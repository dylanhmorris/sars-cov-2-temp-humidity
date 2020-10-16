##################################################
## constants
##################################################

###########################################
## physical constants
###########################################

#' R
#'
#' the ideal gas constant in J per Kelvin mol
#' 
#' @export
R <- 8.31446261815324 # in J per Kelvin mol


###########################################
## unit coversion
###########################################

#' to_Kelvin
#'
#' convert a temperature in Celsius to Kelvin
#'
#' @param temp_C temperature in Celsius
#' @return temperature in Kelvin
#' @export
to_Kelvin <- function(temp_C){

    return(temp_C + 273.15)
}

#' to_Celsius
#'
#' convert a temperature in Kelvin to Celsius
#'
#' @param temp_K temperature in Kelvin
#' @return temperature in Celsius
#' @export
to_Celsius <- function(temp_K){

    return(temp_K - 273.15)
}



###########################################
## empirical equation parameters
###########################################
## parameters for empirical Tang polynomial
a0_tang_NaCl <- 1
a1_tang_NaCl <- (-6.366e-3)
a2_tang_NaCl <- (8.624e-5)
a3_tang_NaCl <- (-1.158e-5)
a4_tang_NaCl <- (1.518e-7)

tang_NaCl_vec <- c(a0_tang_NaCl,
                   a1_tang_NaCl,
                   a2_tang_NaCl,
                   a3_tang_NaCl,
                   a4_tang_NaCl)

## max and min RH for the Tang equation
## to be valid
tang_bounds_NaCl <- c(0.37296, 0.9999999)

## empirical parameters for Pitzer equation
B_0_NaCl <- 0.07722
B_1_NaCl <- 0.25183
C_NaCl <- 0.00106
