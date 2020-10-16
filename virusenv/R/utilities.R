###########################################
## data post-processing
###########################################

#' add_titer_metadata
#'
#' add metadata to a set of MCMC chains
#' ordered by titer_id
#'
#' @param mcmc_draws mcmc chains in spread_draws tibble
#' format from tidybayes
#' @param fitting_data data used to fit the mcmc model
#' as a tibble with a titer id column
#' @param id_column name of the id column, default "titer_id"
#' @return tibble of draws with metadata added
#' @export
add_titer_metadata <- function(mcmc_draws,
                               fitting_data,
                               id_column = "titer_id"){
    
    filt_dat <- dplyr::distinct(fitting_data,
                                fitting_data[id_column],
                                .keep_all = TRUE)
    
    result <- dplyr::inner_join(filt_dat,
                                mcmc_draws,
                                by = id_column)

    if(!dim(result)[1] == dim(mcmc_draws)[1])
        stop(sprintf(paste0("Error: joining data to draws on column %s",
                            "failed to produce set of draws ",
                            "with the same number of rows, ",
                            "suggesting draws were lost; ",
                            "check your data to make sure the ",
                            "join column has all the needed values"),
                     id_column))
                           
    return(result)
}
