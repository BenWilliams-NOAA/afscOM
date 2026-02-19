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
    # return(fy*selex*ret)
    return(sweep(selex*ret, 5, array(fy), "*"))
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
    # return(selex*(1-ret)*dmr)
    return(sweep(selex*(1-ret), 5, array(dmr), "*"))
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

    nfleets <- dim(faa)[length(dim(faa))]

    pred_catches <- array(NA, dim=dim(faa))

    for(f in 1:nfleets) {
        # F-at-age for this fleet
        FAA <- subset_matrix(faa, f, d=5, drop=TRUE)

        # Total Z includes F from ALL fleets
        ZAA_total <- mort
        for(ff in 1:nfleets) {
            ZAA_total <- ZAA_total + subset_matrix(faa, ff, d=5, drop=TRUE)
        }

        # Predicted catch for this fleet (Baranov catch equation)
        pred_catches[,,,,f] <- (FAA / ZAA_total * naa * (1 - exp(-ZAA_total))) * waa
    }
    # return((faa / zaa) * naa * (1 - exp(-zaa)) * waa)
    return(pred_catches)
}
