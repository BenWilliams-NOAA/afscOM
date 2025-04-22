#' Retained Fishing Mortality Rate
#'
#' Calculates the instantaneous fishing-mortality-at-age
#'
#' @param fy the instantaneous fishing mortality rate (F)
#' @param selex selectivity-at-age
#' @param ret retention-at-age
#'
#' @return fishing-mortality-at-age (fy*selex*ret)
#'
#' @export
#'
retained_F <- function(fy, selex, ret){
    return(fy*selex*ret)
}

#' Discarded Fishing Mortality Rate
#'
#' Calculates the instantaneous discard-mortality-at-age
#'
#' @param dmr the instanteous discard mortality rate
#' @param selex selectivity-at-age
#' @param ret retention-at-age
#'
#' @return discard-mortality-at-age (selex*(1-ret)*dmr)
#'
#' @export
#'
discard_F <- function(dmr, selex, ret){
    return(selex*(1-ret)*dmr)
}

#' Catch-at-age
#'
#' Computes catch-at-age given a vector of
#' fishing-mortality-at-age.
#'
#' @param faa instantaneous fishing-mortality-at-age
#' @param naa numbers-at-age vector
#' @param waa weight-at-age vector
#' @param mort instantaneous natural-mortality-at-age
#'
#' @return catch-at-age vector
#'
#' @export
#'
catch_at_age <- function(faa, naa, waa, mort){
    if(all(sapply(list(faa, naa, waa, mort), length) != length(naa))){
        stop("One or more inputs is of different dimensions.")
    }
    zaa <- faa + mort
    return((faa / zaa) * naa * (1 - exp(-zaa)) * waa)
}
