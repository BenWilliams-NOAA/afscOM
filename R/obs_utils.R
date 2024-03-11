#' Simulate Relative Population Numbers Index
#' 
#' Generate a relative population number index based on
#' a survey occurring halfway through the year.
#'
#' @param q catchability coefficient
#' @param naa numbers-at-age vector
#' @param sel selectivity-at-age vector
#' @param zaa total mortality-at-age vector
#'
#' @export simulate_rpn
#'
#' @example
#'
simulate_rpn <- function(q, naa, sel, zaa){
    #total_naa <- apply(naa, 2, sum)
    return(q*sum(naa*exp(-zaa/2)*sel))
}

#' Simulate Relative Population Weight Index
#' 
#' Generate a relative population weight index based on
#' a survey occurring halfway through the year.
#'
#' @param q catchability coefficient
#' @param naa numbers-at-age vector
#' @param waa weight-at-age vector
#' @param sel selectivity-at-age vector
#' @param zaa total mortality-at-age vector
#'
#' @export simulate_rpw
#'
#' @example
#'
simulate_rpw <- function(q, naa, waa, sel, zaa){
    #total_naa <- apply(naa, 2, sum)
    return(q*sum(naa*exp(-zaa/2)*sel*waa))
}

#' Simulate survey age composition data
#' 
#' Generate an age composition vector based on
#' a survey occurring halfway through the year,
#' and for which ageing error may exist.
#'
#' @param naa numbers-at-age vector
#' @param sel selectivity-at-age vector
#' @param age_err aging-error matrix
#'
#' @export simulate_ac
#'
#' @example
#'
simulate_ac <- function(naa, sel, age_err=NA){
    eac <- naa*sel
    eac <- apply(eac, c(1, 2), sum)
    if(!all(is.na(age_err))){
        eac <- eac %*% age_err
    }
    return(eac/sum(eac))
}
