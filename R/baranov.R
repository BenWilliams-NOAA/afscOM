#' Baranov Catch Equation
#'
#' An implementation of the Baranov catch equations allowing for
#' fishery retention and discarding practices.
#'
#' @param fy an instantaneous fishing mortality rate to apply to the population
#' @param naa numbers-at-age in the population
#' @param waa weight-at-age vector
#' @param mort instantaneous natural mortality (M) -at-age vector
#' @param selex selectivity-at-age vector
#' @param ret retention-at-age vector (optional)
#' @param dmr instantaneous discard mortality rate -at-age vector (optional)
#'
#' @return the amount of catch resulting from fishing the population at `fy`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' naa <- c(100, 50, 25, 10)
#' waa <- c(1, 2, 3, 4)
#' mort <- c(0.2, 0.2, 0.2, 0.2)
#' selex <- c(0.0, 0.5, 1.0, 1.0)
#'
#' # catch should be 14.35
#' catch <- baranov(0.1, naa, waa, mort, selex)
#'}
baranov <- function(fy, naa, waa, mort, selex, ret=NA, dmr=NA) {
    if(any(is.na(ret)) | any(is.na(dmr))){
        warning("No discard mortality rate OR retention-at-age provided.\nAssuming full retention.")
        ret <- rep(1, length(selex))
        dmr <- rep(0, length(ret))
    }

    if(all(sapply(list(naa, waa, mort, selex, ret, dmr), length) != length(naa))){
        stop("One or more inputs is of different dimensions.")
    }

    faa <- retained_F(fy, selex, ret) + discard_F(dmr, selex, ret)
    zaa <- faa + mort
    caa <- naa*faa*(1-exp(-zaa))/zaa
    catch <- sum(caa*waa)
    return(catch)
}
