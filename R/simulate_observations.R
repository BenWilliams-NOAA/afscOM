#' Simulate Observations from Fisheries and Surveys
#'
#' Simulate observations of a population from common Alaska surveys and 
#' from landed catch data.
#'
#' @param naa numbers-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])
#' @param waa weight-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])
#' @param selex selectivity-at-age array subset to 1 year (dim [1, nages, nsexes, nregions, nfleets])
#' @param faa fishing mortality-at-age array subset to 1 year (dim [1, nages, nsexes, nregions, nfleets])
#' @param zaa total mortality-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])
#' @param obs_pars list of observation process parameters
#' @param age_error an optional ageing error matrix
#'
#' @export simulate_observations
#'
simulate_observations <- function(naa, waa, selex, faa, zaa, obs_pars, age_error=NA){

    ll_selex <- subset_matrix(selex, r=1, d=5, drop=TRUE)
    tw_selex <- subset_matrix(selex, r=2, d=5, drop=TRUE)

    fxfish_faa <- subset_matrix(faa, r=1, d=5, drop=TRUE)
    twfish_faa <- subset_matrix(faa, r=2, d=5, drop=TRUE)

    model_params <- get_model_dimensions(selex)

    # Simulate Domestic LL Survey RPN (log-normal)
    # --------------------------------------------
    ll_rpn_pred <- simulate_rpn(obs_pars$surv_ll$q, naa, ll_selex, zaa)
    ll_rpn_obs  <- simulate_lognormal_obs(ll_rpn_pred, obs_pars$surv_ll$rpn_cv)

    ll_rpw_pred <- simulate_rpw(obs_pars$surv_ll$q, naa, waa, ll_selex, zaa)
    ll_rpw_obs  <- simulate_lognormal_obs(ll_rpw_pred, obs_pars$surv_ll$rpw_cv)

    # Simulate Trawl Survey RPW (log-normal)
    # --------------------------------------
    tw_rpw_pred <- simulate_rpw(obs_pars$surv_tw$q, naa, waa, tw_selex, zaa)
    tw_rpw_obs  <- simulate_lognormal_obs(tw_rpw_pred, obs_pars$surv_tw$rpw_cv)

    # Simulate Domestic LL Survey Age Compositions (multinomial)
    # ----------------------------------------------------------
    ll_ac_pred <- simulate_ac(naa, ll_selex, aggregate_sex = FALSE)
    # Perform multinomial draw for each sex
    ll_ac_obs <- simulate_multinomial_obs(ll_ac_pred, obs_pars$surv_ll$ac_samps, as_integers = obs_pars$surv_ll$as_integers)
    ll_ac_obs <- array(ll_ac_obs, dim=c(model_params$nyears, model_params$nages, length(ll_ac_obs)/model_params$nages, model_params$nregions), dimnames=dimnames(naa))

    # Simulate Trawl Survey Age Compositions (multinomial)
    # ----------------------------------------------------
    tw_ac_pred <- simulate_ac(naa, tw_selex, aggregate_sex = FALSE)
    tw_ac_obs <- simulate_multinomial_obs(tw_ac_pred, obs_pars$surv_tw$ac_samps, as_integers = obs_pars$surv_tw$as_integers)
    tw_ac_obs <- array(tw_ac_obs, dim=c(model_params$nyears, model_params$nages, length(tw_ac_obs)/model_params$nages, model_params$nregions), dimnames=dimnames(naa))


    # Simulate Domestic LL Fishery Age Compositions (multinomial)
    # -----------------------------------------------------------
    fxfish_caa_pred <- simulate_caa(naa, fxfish_faa, zaa)
    fxfish_caa_obs  <- simulate_multinomial_obs(fxfish_caa_pred, obs_pars$fish_fx$ac_samps, as_integers = obs_pars$fish_fx$as_integers, age_err=NA)
    fxfish_caa_obs  <- array(fxfish_caa_obs, dim=c(model_params$nyears, model_params$nages, length(fxfish_caa_obs)/model_params$nages, model_params$nregions), dimnames=dimnames(naa))

    # Simulate Domestic Trawl Fishery Age Compositions (multinomial)
    # -----------------------------------------------------------
    twfish_caa_pred <- simulate_caa(naa, twfish_faa, zaa)
    twfish_caa_obs  <- simulate_multinomial_obs(twfish_caa_pred, obs_pars$fish_tw$ac_samps, as_integers = obs_pars$fish_tw$as_integers, age_err=NA)
    twfish_caa_obs  <- array(twfish_caa_obs, dim=c(model_params$nyears, model_params$nages, length(twfish_caa_obs)/model_params$nages, model_params$nregions), dimnames=dimnames(naa))

    # Simulate Trawl Fishery Age Compositions (mulitnomial)
    preds <- listN(ll_rpn_pred, ll_rpw_pred, tw_rpw_pred, ll_ac_pred, tw_ac_pred, fxfish_caa_pred, twfish_caa_pred)
    obs   <- listN(ll_rpn_obs, ll_rpw_obs, tw_rpw_obs, ll_ac_obs, tw_ac_obs, fxfish_caa_obs, twfish_caa_obs)

    return(list(
        preds=preds,
        obs=obs
    ))

}
