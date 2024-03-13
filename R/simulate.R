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
#' @export sim_rpn
#'
#'
sim_rpn <- function(q, naa, sel, zaa){
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
#' @export sim_rpw
#'
#'
sim_rpw <- function(q, naa, waa, sel, zaa){
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
#' @export sim_ac
#'
#'
sim_ac <- function(naa, sel, aggregate_sex=FALSE){
    eac <- naa*sel
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
#' @export sim_caa
#'
#'
sim_caa <- function(naa, faa, zaa, aggregate_sex=FALSE){

    caa <- naa*faa*(1-exp(-zaa))/zaa
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
#' @export sim_lognormal_obs
#'
#' @example sim_lognormal_obs(10, 0.20)
#'
sim_lognormal_obs <- function(pred, cv){
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
#' @export sim_multinomial_obs
#'
#' @example sim_multinomial_obs(c(0.25, 0.50, 0.25), 100)
#'
sim_multinomial_obs <- function(pred, samp_size, age_err=NA){
    multi <- apply(pred, 3, function(x){
                tmp <- rmultinom(1, samp_size, prob=x)
                tmp <- (tmp[,1]/sum(tmp[,1]))
                # Include aging error if available
                if(!all(is.na(age_err))){
                    tmp <- tmp %*% age_err
                }
                return(tmp)
            })
    return(multi)
}

sim_pop <- function(prev_naa, faa, recruitment, dem_params, options){

    model_params <- get_model_dimensions(dem_params$sel)

    faa <- array(apply(faa, c(2, 3), sum), dim=c(1, model_params$nages, model_params$nsexes, 1)) # Sum across fleets to get total FAA
    Zaa <- faa + dem_params$mort

    naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))

    # Recruitment goes here
    naa[1,1,,1] <- recruitment

    # Simulate mortality here
    for(a in 2:(model_params$nages-1)){
        naa[1, a, ,1] <- prev_naa[1, a-1,,1]*exp(-Zaa[,a-1,,])
    }
    naa[1, model_params$nages, ,1] <- prev_naa[1, model_params$nages-1,,1]*exp(-Zaa[,model_params$nages-1,,])+prev_naa[1, model_params$nages,,1]*exp(-Zaa[,model_params$nages,,])

    return(listN(naa))

}

#' Simulate Observations from Fisheries and Surveys
#'
#' Description
#'
#' @param $$
#'
#' @export
#'
sim_obs <- function(naa, waa, selex, faa, zaa, obs_pars, age_error=NA){

    ll_selex <- subset_matrix(selex, r=1, d=5, drop=TRUE)
    tw_selex <- subset_matrix(selex, r=2, d=5, drop=TRUE)

    fxfish_faa <- subset_matrix(faa, r=1, d=5, drop=TRUE)
    twfish_faa <- subset_matrix(faa, r=2, d=5, drop=TRUE)

    model_params <- get_model_dimensions(selex)

    # Simulate Domestic LL Survey RPN (log-normal)
    # --------------------------------------------
    ll_rpn_pred <- sim_rpn(obs_pars$surv_ll$q, naa, ll_selex, zaa)
    ll_rpn_obs  <- sim_lognormal_obs(ll_rpn_pred, obs_pars$surv_ll$rpn_cv)

    ll_rpw_pred <- sim_rpw(obs_pars$surv_ll$q, naa, waa, ll_selex, zaa)
    ll_rpw_obs  <- sim_lognormal_obs(ll_rpw_pred, obs_pars$surv_ll$rpw_cv)

    # Simulate Trawl Survey RPW (log-normal)
    # --------------------------------------
    tw_rpw_pred <- sim_rpw(obs_pars$surv_tw$q, naa, waa, tw_selex, zaa)
    tw_rpw_obs  <- sim_lognormal_obs(tw_rpw_pred, obs_pars$surv_tw$rpw_cv)

    # Simulate Domestic LL Survey Age Compositions (multinomial)
    # ----------------------------------------------------------
    ll_ac_pred <- sim_ac(naa, ll_selex, aggregate_sex = FALSE)
    # Perform multinomial draw for each sex
    ll_ac_obs <- sim_multinomial_obs(ll_ac_pred, obs_pars$surv_ll$ac_samps)
    ll_ac_obs <- array(ll_ac_obs, dim=c(model_params$nyears, model_params$nages, 
                        length(ll_ac_obs)/model_params$nages, model_params$nregions), 
                        dimnames=dimnames(naa))

    # Simulate Trawl Survey Age Compositions (multinomial)
    # ----------------------------------------------------

    # Simulate Domestic LL Fishery Age Compositions (multinomial)
    # -----------------------------------------------------------
    fxfish_caa_pred <- sim_caa(naa, fxfish_faa, zaa)
    fxfish_caa_obs  <- sim_multinomial_obs(fxfish_caa_pred, obs_pars$fish_fx$ac_samps, age_err=NA)
    fxfish_caa_obs <- array(fxfish_caa_obs, dim=c(model_params$nyears, model_params$nages, 
                            length(fxfish_caa_obs)/model_params$nages, model_params$nregions), 
                            dimnames=dimnames(naa))

    # Simulate Trawl Fishery Age Compositions (mulitnomial)
    preds <- listN(ll_rpn_pred, ll_rpw_pred, tw_rpw_pred, ll_ac_pred, fxfish_caa_pred)
    obs   <- listN(ll_rpn_obs, ll_rpw_obs, tw_rpw_obs, ll_ac_obs, fxfish_caa_obs)

    return(list(
        preds=preds,
        obs=obs
    ))

}

    #' Simulate catches in one region during on year
#'
#' Description
#'
#' @param TAC the total catch to be removed across all fleets
#' @param fleet_props the proportion of the TAC that should be taken by each fleet
#' @param dem_params demographic parameter matrices for population
#' @param naa the current numbers-at-age in the population
#' @param option list of model options
#'
#' @export
#'
sim_catch <- function(TAC, fleet_props, dem_params, naa, options){

    model_params <- get_model_dimensions(dem_params$sel)

    # if not provided, assume equal fleet split
    if(!("fleet_apportionment" %in% names(options))){
        options$fleet_apportionment <- rep(1/model_params$nfleets, model_params$nfleets)
    }

    land_caa_tmp    <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    disc_caa_tmp    <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    caa_tmp         <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    faa_tmp         <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))

    for(f in 1:model_params$nfleets){

        # Skip if the fleet in question is a survey fleet
        # if("fleet_is_survey" %in% model_options & model_options$fleet_is_survey[f]){
        #     next;
        # }

        tac <- TAC*fleet_props[f]

        dp_f <- rlang::duplicate(dem_params)
        dp_f$sel <- subset_matrix(dp_f$sel, f, d=5, drop=TRUE)
        dp_f$ret <- subset_matrix(dp_f$ret, f, d=5, drop=TRUE)
        dp_f$dmr <- subset_matrix(dp_f$dmr, f, d=5, drop=TRUE)

        F_f <- findF_bisection(
            f_guess = 0.05,
            naa     = naa,
            waa     = dp_f$waa,
            mort    = dp_f$mort,
            selex   = dp_f$sel,
            ret     = dp_f$ret,
            dmr     = dp_f$dmr,
            prov_catch = tac
        )

        # ret_faa <- retained_F(F_f, dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        # disc_faa <- discard_F(dem_params$dmr[,,,,f], dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        ret_faa <- retained_F(F_f, dp_f$sel, dp_f$ret)
        disc_faa <- discard_F(dp_f$dmr, dp_f$sel, dp_f$ret)
        faa_tmp[,,,,f] <- ret_faa + disc_faa

        land_caa_tmp[,,,,f] <- catch_at_age(ret_faa, naa, dp_f$waa, dp_f$mort)
        disc_caa_tmp[,,,,f] <- catch_at_age(disc_faa, naa, dp_f$waa, dp_f$mort)
        caa_tmp[,,,,f] <- land_caa_tmp[,,,,f] + disc_caa_tmp[,,,,f]
    }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, faa_tmp))

}