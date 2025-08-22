#' Compute Spawning Stock Biomass
#'
#' Compute SSB from a matrix of NAA and a set of demographic
#' parameters.
#'
#' @param naa numbers-at-age matrix [nyears, nages, nsexes, nregions]
#' @param dem_params demagraphic parameters list
#'
#' @export compute_ssb
#'
compute_ssb <- function(naa, dem_params){
    return(
        apply(naa[1:nrow(dem_params$waa),,1,,drop=FALSE]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], c(1, 4), sum)
    )
}

#' Compute Total Stock Biomass
#'
#' Compute biomass from a matrix of NAA and a set of demographic
#' parameters. This computes the sum of male and female biomass
#' across all age classes. For mature female biomass, see
#' `compute_ssb`.
#'
#' @param naa numbers-at-age matrix [nyears, nages, nsexes, nregions]
#' @param dem_params demagraphic parameters list
#'
#' @export compute_bio
#'
compute_bio <- function(naa, dem_params){
    return(
        apply(naa[1:nrow(dem_params$waa),,,,drop=FALSE]*dem_params$waa[,,,,drop=FALSE], c(1, 4), sum)
    )
}

#' Compute Total Catch
#'
#' Compute total catch volume (in mt) from a matrix of CAA.
#'
#' @param caa catch-at-age matrix [nyears, nages, nsexes, nregions, nfleets]
#'
#' @export compute_total_catch
#'
compute_total_catch <- function(caa){
    return(
        apply(caa, c(1, 4), sum)
    )
}

#' Compute Catch by Fishing Fleet
#'
#' Compute total catch volume (in mt) for each fishing
#' fleet from a matrix of CAA.
#'
#' @param caa catch-at-age matrix [nyears, nages, nsexes, nregions, nfleets]
#'
#' @export compute_fleet_catch
#'
compute_fleet_catch <- function(caa){
    return(
        apply(caa, c(1, 4, 5), sum)
    )
}

#' Compute Total Fishing Mortality
#'
#' Compute total fishing mortality from a matrix of FAA.
#' Note that this is an approximation of total fishing mortality,
#' and assumes that selectivity = 1 for some age class.
#'
#' @param faa fishing-mortality-at-age matrix [nyears, nages, nsexes, nregions, nfleets]
#'
#' @export compute_total_f
#'
compute_total_f <- function(faa){
    return(
        apply(apply(faa, c(1, 4, 5), \(x) max(x)), c(1, 2), sum)
    )
}
