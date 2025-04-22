#' Find F given a desired level of catch
#'
#' Uses the Baranov catch equation to compute the instantaneous
#' fishing mortality rate (F) that would yield a specific level of
#' catch_
#'
#' @param f_guess initial guess for F
#' @param naa numbers-at-age
#' @param waa weight-at-age
#' @param mort natural-mortality-at-age
#' @param selex selectivity-at-age
#' @param ret retention-at-age (optional)
#' @param dmr instantaneous discard mortality rate (optional)
#' @param prov_catch desired catch level
#'
#' @return the fishing mortality rate that yields the desired catch level
#'
#' @export
#'
#' @examples
#' \dontrun{
#' find_F(f_guess=0.1, naa=naa, waa=waa, mort=mort, selex=selex, prov_catch=0.12)
#' }
#'
find_F <- function(f_guess, naa, waa, mort, selex, ret=NA, dmr=NA, prov_catch){

    # Quick switch so as to not run the optimization routine
    # if no catch is desired.
    if(prov_catch == 0){
        return(0.0)
    }

    baranov_min <- function(ln_fy, naa, waa, mort, selex, ret, dmr, prov_catch) {
        fy = exp(ln_fy)
        catch <- baranov(fy, naa, waa, mort, selex, ret, dmr)
        min <- abs(prov_catch - catch)
        return(min)
    }

    # use downhill simplex to get it (no need for derivatives which is great for minimizing abs vals)
    # NOTE: ideally, we eliminate the dependency on bbmle::mle2 and replace this with a bisection
    #       algorithm.
    c_to_f = bbmle::mle2(baranov_min,
        start = list(ln_fy=log(f_guess)),
        vecpar = TRUE,
        data = list(
            naa = naa,
            waa = waa,
            mort = mort,
            selex = selex,
            ret = ret,
            dmr = dmr,
            prov_catch = prov_catch
        ),
        method="Nelder-Mead",
        optimizer="nlminb",
        control=list(eval_max=5e6, iter_max=5e6)
    )
    return(exp(c_to_f@coef))
}

#' Find F given a desired level of catch
#'
#' Uses the Baranov catch equation to compute the instantenous
#' fishing mortality rate (F) that would yield a specific level of
#' catch using bisection
#'
#' @param f_guess initial guess for F
#' @param naa numbers-at-age
#' @param waa weight-at-age
#' @param mort natural-mortality-at-age
#' @param selex selectivity-at-age
#' @param ret retention-at-age (optional)
#' @param dmr instantaneous discard mortality rate (optional)
#' @param prov_catch desired catch level
#'
#' @return the fishing mortality rate that yields the desired catch level
#'
#' @export findF_bisection
#'
#' @examples
#' \dontrun{
#' findF_bisection(f_guess=0.1, naa=naa, waa=waa, mort=mort, selex=selex, prov_catch=12)
#' }
#'
findF_bisection <- function(f_guess, naa, waa, mort, selex, ret=NA, dmr=NA, prov_catch){

    range <- vector(length=2)
    range[1] <- 0
    range[2] <- 1
    n.iter <- 20
    i <- 1
    for(i in 1:n.iter) {
      midpoint <- mean(range)
      pred_catch <- baranov(midpoint, naa, waa, mort, selex, ret, dmr)
      if(pred_catch < prov_catch) {
        range[1] <- midpoint
        range[2] <- range[2]
      }else {
        range[1] <- range[1]
        range[2] <- midpoint
      }
    }
    Fmort <- midpoint
    return(Fmort)
}
