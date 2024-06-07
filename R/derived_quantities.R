#' Compute Spawning Stock Biomass
#' 
#' Compute SSB from a matrix of NAA and a set of demographic
#' parameters.
#'
#' @param naa numbers-at-age matrix
#' @param dem_params demagraphic parameters list
#'
#' @export compute_ssb
#'
#' @example
#'
compute_ssb <- function(naa, dem_params){
    return(
        apply(naa[1:nrow(dem_params$waa),,1,,drop=FALSE]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], c(1, 4), sum)
    )
}

compute_bio <- function(naa, dem_params){
    return(
        apply(naa[1:nrow(dem_params$waa),,,,drop=FALSE]*dem_params$waa[,,,,drop=FALSE], c(1, 4), sum)
    )
}

compute_total_catch <- function(caa){
    return(
        apply(caa, c(1, 4), sum)
    )
}

compute_fleet_catch <- function(caa){
    return(
        apply(caa, c(1, 4, 5), sum)
    )
}

compute_total_f <- function(faa){
    return(
        apply(apply(faa, c(1, 4, 5), \(x) max(x)), c(1, 2), sum)
    )
}
