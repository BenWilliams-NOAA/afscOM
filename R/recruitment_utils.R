#' Beverton-Holt Stock Recruit Relationship
#' 
#' Classic beverton-holt SRR parameterized using steepness,
#' unfished recruitment (R0), and unfished spawning biomass
#' (S0). 
#'
#' @param naa numbers-at-age array (dimensions [1, nages, nsexes, nregions])
#' @param dem_params demographic parameters list susbet to a 
#' single year (dimensions [1, nages, nsexes, nregions, nfleets])
#'
#' @export beverton_holt
#'
#' @example
#'
beverton_holt <- function(naa, dem_params, h, R0, S0, sigR, seed){
    set.seed(seed)
    ssb <- compute_ssb(naa, dem_params)[1,1]
    bh <- (4*R0*h*ssb)/((1-h)*R0*(S0/R0) + (5*h - 1)*ssb)
    rec <- rlnorm(1, meanlog = log(bh)-sigR*sigR/2, sdlog=sigR)
    return(rec)
}

