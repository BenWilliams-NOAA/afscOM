#' Resample from historical recruitment timeseries
#'
#' Draw `n` random samples, with replacement, from the
#' historical recruitment timeseries.
#'
#' @param n number of samples to draw
#' @param hist_rec historical recruitment vector
#' @param seed a random seed (optional)
#'
#' @return `n` random samples from the historical
#' recruitment timeseries.
#'
#' @export
#'
#' @examples
#' set.seed(1120)
#' n <- 100
#' hist.rec <- c(100, 200, 50, 100, 120)
#' resample.recruitment(n, hist.rec)
#'
resample_recruitment <- function(n, hist_rec, seed=NA){
    if(!is.na(seed)){
        set.seed(seed)
    }
    return(sample(hist_rec, n, replace=TRUE))
}

mean_recruitment <- function(n, log.mu, sd, seed=NA){
    if(!is.na(seed)){
        set.seed(seed)
    }
    return(exp(log.mu + stats::rnorm(n, 0, sd) - sd^2/2))
}

#' Generate log-recruitment from two defined regimes of a specified length
#'
#' Description
#'
#' @param n number of recruitment events
#' @param mu_rec vector of mean log-recruitment levels, one for each regime
#' @param sd_rec vector of recruitment variation, one for each regime
#' @param dur duration of regimes
#' @param start which regime to start in
#' @param seed a random seed (optional)
#'
#' @return `n` log-recruitment levels from two different regimes of
#' length `dur`.
#' @export
#'

cyclic_recruitment <- function(n, mu_rec, sd_rec, dur, start=1, seed=NA){

    if(!is.na(seed)){
        set.seed(seed)
    }

    curr_regime <- start
    duration <- dur[curr_regime]
    recs <- rep(NA, n)
    for(i in 0:n){
        recs[i] <- exp(mu_rec[curr_regime+1] + rnorm(1, 0, sd_rec[curr_regime+1]) - (sd_rec[curr_regime+1]^2)/2)
        if(i %% duration == 0){
            curr_regime <- as.numeric(!curr_regime)
            duration <- dur[curr_regime+1]
        }
    }
    return(recs)

}

#' Beverton-Holt Stock Recruit Relationship
#'
#' Classic Beverton-Holt SRR parameterized using steepness,
#' unfished recruitment (R0), and unfished spawning biomass
#' (S0).
#'
#' @param naa numbers-at-age array (dimensions [1, nages, nsexes, nregions])
#' @param dem_params demographic parameters list susbet to a
#' single year (dimensions [1, nages, nsexes, nregions, nfleets])
#'
#' @export
#'
#' @examples
#' \dontrun{
#' beverton_holt(naa, dem_params, h, R0, S0, sigR, seed)
#' }
#'
#'
beverton_holt <- function(naa, dem_params, h, R0, S0, sigR, seed){
  set.seed(seed)
  ssb <- compute_ssb(naa, dem_params)[1,1]
  bh <- (4*R0*h*ssb)/((1-h)*R0*(S0/R0) + (5*h - 1)*ssb)
  rec <- stats::rlnorm(1, meanlog = log(bh)-sigR*sigR/2, sdlog=sigR)
  return(rec)
}

