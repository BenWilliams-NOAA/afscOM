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
#' @param aggregate_sex whether to return sex-aggregated comps or not (default FALSE
#'
#' @export simulate_ac
#'
#'
simulate_ac <- function(naa, sel, aggregate_sex=FALSE){
    eac <- naa*sel
    if(all(eac == 0)){
        return(array(0, dim=dim(naa), dimnames=dimnames(naa)))
    }

    if(aggregate_sex){
        eac <- array(apply(eac, c(1, 2), sum), dim=c(1, dim(naa)[2], 1))
    }
    
    # if(!all(is.na(age_err))){
    #     eac <- eac %*% age_err
    # }

    if(aggregate_sex){ 
        out_dims <- c(1, dim(naa)[2], 1) 
    }else{ 
        out_dims <- dim(naa)
    }
    std_comp <- array(apply(eac, c(3), \(x) x/sum(x)), dim=out_dims)

    return(std_comp)
}

#' Simulate catch-at-age composition data
#'
#' Generate an age composition vector based on
#' fishing induced mortality.
#'
#' @param naa numbers-at-age vector
#' @param faa fishing-mortality-at-age vector
#' @param zaa total-mortality-at-age vector
#' @param aggregate_sex whether to return sex-aggregated comps or not (default FALSE)
#'
#' @export simulate_caa
#'
#'
simulate_caa <- function(naa, faa, zaa, aggregate_sex=FALSE){

    caa <- naa*faa*(1-exp(-zaa))/zaa

    # Short circuit, return 0s is no catch
    if(all(caa == 0)){
        return(array(0, dim=dim(caa), dimnames=dimnames(caa)))
    }

    caa_prop <- array(apply(caa, c(3), \(x) x/rowSums(x)), dim=dim(naa), dimnames=dimnames(naa))
    eac <- caa_prop
    if(aggregate_sex){
        eac <- array(apply(caa_prop, c(1, 2), sum)/2, dim=c(1, dim(naa)[2], 1))
    }

    # if(!all(is.na(age_err))){
    #     eac <- eac %*% age_err
    # }

    if(aggregate_sex){ 
        out_dims <- c(1, dim(naa)[2], 1) 
    }else{ 
        out_dims <- dim(naa)
    }
    std_comp <- array(apply(eac, c(3), \(x) x/sum(x)), dim=out_dims)

    return(std_comp)
}

#' Simulate observations from a lognormal distribution
#' #'
#' A wrapper function around `rlnorm` that generates a
#' single random observation from a lognormal distribution
#' centered on a predicted value (`pred`) and given a 
#' level of error (`cv`).
#'
#' @param pred the predicted value of an observation
#' @param cv the real-space coefficient of variation about
#' the predicted observation value
#'
#' @export simulate_lognormal_obs
#'
#' @example simulate_lognormal_obs(10, 0.20)
#'
simulate_lognormal_obs <- function(pred, cv){
    sds <- sqrt(log(cv^2 + 1))
    return(rlnorm(1, meanlog=log(pred)-(sds^2)/2, sdlog = sds))
}

#' Simulate observations from a multinomial distribution
#' #'
#' A wrapper function around `rmultinom` that generates a
#' single vector of observations from a multinomial distribution
#' where each of the `K` multinomial classes has probability
#' `pred`. This function will automatically handle generating
#' multinomial draws across multiple sexes. Aging error can be
#' optionally applied. 
#'
#' @param pred the probability associated with each class
#' @param samp_size multinomical sample size
#' @param age_err whether to automaticaly apply ageing error
#'
#' @export simulate_multinomial_obs
#'
#' @example simulate_multinomial_obs(c(0.25, 0.50, 0.25), 100)
#'
simulate_multinomial_obs <- function(pred, samp_size, age_err=NA){
    multi <- array(0, dim=dim(pred), dimnames=dimnames(pred))
    if(!all(pred == 0)){
        multi <- apply(pred, 3, function(x){
                tmp <- rmultinom(1, samp_size, prob=x)
                tmp <- (tmp[,1]/sum(tmp[,1]))
                # Include aging error if available
                if(!all(is.na(age_err))){
                    tmp <- tmp %*% age_err
                }
                return(tmp)
            })
    }
    return(multi)
}
