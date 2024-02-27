#' mu_to_F
#' #'
#' Converts between harvest rate (the proportion
#' of the population harvest; \{mu}) and instantenous
#' fishing mortality rate (F)_
#' 
#' @param mu, a harvest rate
#'
#' @return F, the corresponding instantenous fishing mortality rate
#'
#' @export mu_to_F
#'
mu_to_F <- function(mu){
    return(-log(1-mu))
}

#' F_to_mu
#' 
#' Converts between instantaneous fishing mortality
#' rate (F) and harvest rate (the proportion of the
#' population harvest; \{mu})_
#'
#' @param fy, an instantaneous fishing mortality rate
#'
#' @return mu, the corresponding harvest rate
#'
#' @export F_to_mu
#'
F_to_mu <- function(fy){
    return(1-exp(-fy))
}

#' Set variable names as names of list
#' 
#' Automatically use variable names as list object names when
#' a sequence of variables is provided to a list constructor_
#'
#' @param ... R variables containung values to be put in a list
#'
#' @export listN
#'
#' @example 
#' a = b = d = e = 1
#' listN(a, b, d, e)
#'
listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

#' Subset demographic parameter matrices by first dimensions
#' 
#' Access an arbitrary row across all demographic parameter
#' matrices of indeterminate dimension
#'
#' @param dem_params a list of arrays with more than 2 dimensions
#' @param r the index along the first dimension to access
#' 
#' @return a list with the same elements as dem_params but containing
#' a single row from each list element
#'
#' @export subset_dem_params
#'
#' @example 
#'
subset_dem_params <- function(dem_params, r, d=1){
    tmp <- rlang::duplicate(dem_params)
    ps <- names(tmp)
    for(n in ps){
        ndims <- length(dim(dem_params[[n]]))
        idxs <- c(as.list(rep(TRUE, d-1)), list(r), as.list(rep(TRUE, ndims-d)))
        tmp[[n]] <- do.call('[', c(list(dem_params[[n]]), idxs))
    }
    return(tmp)
}