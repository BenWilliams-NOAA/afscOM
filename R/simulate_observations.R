#' Simulate Observations from Fisheries and Surveys
#' 
#' Description
#'
#' @param $$
#'
#' @export
#'
#' @example
#'
simulate_observations <- function(naa, waa, selex, faa, zaa, obs_pars, age_error=NA){

    ll_selex <- subset_matrix(selex, r=1, d=5, drop=TRUE)
    tw_selex <- subset_matrix(selex, r=2, d=5, drop=TRUE)

    fxfish_faa <- subset_matrix(faa, r=1, d=5, drop=TRUE)
    twfish_faa <- subset_matrix(faa, r=2, d=5, drop=TRUE)

    # Simulate Domestic LL Survey RPN (log-normal)
    ll_rpn_pred <- simulate_rpn(obs_pars$surv_ll_q, naa, ll_selex, zaa)
    # ll_rpn_obs  <- rlnorm(1, meanlog = log(ll_rpn_pred), sdlog=obs_pars$ll_rpn_se)

    ll_rpw_pred <- simulate_rpw(obs_pars$surv_ll_q, naa, waa, ll_selex, zaa)
    # ll_rpw_obs  <- rlnorm(1, meanlog = log(ll_rpw_pred), sdlog=obs_pars$ll_rpw_se)

    # Simulate Trawl Survey RPW (log-normal)
    tw_rpw_pred <- simulate_rpw(obs_pars$surv_tw_q, naa, waa, tw_selex, zaa)
    # tw_rpw_obs  <- rlnorm(1, meanlog = log(tw_rpw_pred), sdlog=obs_pars$tw_rpw_se)

    # Simulate Domestic LL Survey Age Compositions (multinomial)
    ll_ac_pred <- simulate_ac(naa, ll_selex, age_err=age_error)
    # ll_ac_obs  <- rmultinom(1, obs_pars$ll_ac_n, prob=ll_ac_pred)
    # ll_ac_obs  <- ll_ac_obs/sum(ll_ac_obs)

    # Simulate Trawl Survey Age Compositions (multinomial)

    # Simulate Domestic LL Fishery Age Compositions (multinomial)
    fxfish_caa_pred <- simulate_caa(naa, fxfish_faa, zaa, age_err = age_error)
    # fxfish_caa_obs  <- rmultinom(1, obs_pars$fxfish_caa_n, prob=fxfish_caa_pred)
    # fxfish_caa_obs  <- lfxfish_caa_obs/sum(fxfish_caa_obs)

    # Simulate Trawl Fishery Age Compositions (mulitnomial)
    preds <- listN(ll_rpn_pred, ll_rpw_pred, tw_rpw_pred, ll_ac_pred, fxfish_caa_pred)
    # obs   <- listN(ll_rpn_obs, ll_rpw_obs, tw_rpw_obs, ll_ac_obs, fxfish_caa_obs)

    return(list(
        preds=preds,
        obs=list()
    ))

}
