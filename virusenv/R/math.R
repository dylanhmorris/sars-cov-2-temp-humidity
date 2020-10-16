#' solve_quartic
#'
#' analytically sole a quartic equation with the given coefficients
#' 
#' @param a_vec vector of coefficients in order a0, a1, a2, a3, a4 \cr
#' where the quartic is a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0
#' @return the four roots
#' @export
solve_quartic <- function(a_vec){
    ## do this explicitly because
    ## R is 1-indexed and we don't
    ## want to screw up!
    e <- as.complex(a_vec[1])
    d <- as.complex(a_vec[2])
    c <- as.complex(a_vec[3])
    b <- as.complex(a_vec[4])
    a <- as.complex(a_vec[5])
    

    p <- (8 * a * c - 3 * b^2) / (8 * a^2)
    q <- (b^3 - 4 * a * b * c + 8 * a^2 * d) / (8 * a^3)

    Delta0 <- c^2 - 3 * b * d + 12 * a * e
    Delta1 <- (2 * c^3 -
               9 * b * c * d +
               27 * b^2 * e +
               27 * a * d^2 -
               72 * a * c * e)

    Q <- ((Delta1 + sqrt(Delta1^2 - 4 * Delta0^3)) / 2)^(1/3)
    
    S <- 0.5 * sqrt((-2 * p / 3) +
                    (1 / (3 * a)) *
                    (Q + Delta0 / Q))
    
    term1 <- -b / (4 * a)
    term3_plus <- 0.5 * sqrt(-4 * S^2 - 2 * p + q / S)
    term3_minus <- 0.5 * sqrt(-4 * S^2 - 2 * p - q / S)
    x1 <- term1 - S + term3_plus
    x2 <- term1 - S - term3_plus
    x3 <- term1 + S + term3_minus
    x4 <- term1 + S - term3_minus
    
    return(c(x1, x2, x3, x4))
}


#' solve_tang_quartic
#'
#' analytically sole a Tang quartic equation with
#' the given coefficients using the fact that the
#' given root will always be the smallest real
#' root on the interval of interest 
#' 
#' @param a_vec vector of coefficients in order a0, a1, a2, a3, a4 \cr
#' where the quartic is a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0
#' @return the four roots
#' @export
solve_tang_quartic <- function(a_vec){
    ## do this explicitly because
    ## R is 1-indexed and we don't
    ## want to screw up!
    e <- a_vec[1]
    d <- a_vec[2]
    c <- a_vec[3]
    b <- a_vec[4]
    a <- a_vec[5]
    

    p <- (8 * a * c - 3 * b^2) / (8 * a^2)
    q <- (b^3 - 4 * a * b * c + 8 * a^2 * d) / (8 * a^3)

    Delta0 <- c^2 - 3 * b * d + 12 * a * e
    Delta1 <- (2 * c^3 -
               9 * b * c * d +
               27 * b^2 * e +
               27 * a * d^2 -
               72 * a * c * e)

    Q <- ((Delta1 + sqrt(Delta1^2 - 4 * Delta0^3)) / 2)^(1/3)
    
    S <- 0.5 * sqrt((-2 * p / 3) +
                    (1 / (3 * a)) *
                    (Q + Delta0 / Q))
    
    term1 <- -b / (4 * a)
    term3_minus <- 0.5 * sqrt(-4 * S^2 - 2 * p - q / S)
    x4 <- term1 + S - term3_minus
    
    return(x4)
}
