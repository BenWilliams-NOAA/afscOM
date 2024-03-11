#' Project population forward by 1 year
#'
#' Given a TAC, the previous population age structure, demographic information
#' (e.g. natural mortality, selectivity, retention, etc.), and recruitment,
#' determine the numbers-at-age present in the population for the following year.
#'
#' @param TAC the total catch to be removed from the population
#' @param fleet.props the proportion of the TAC allocated to each fleet
#' @param dem_params list of demographic parameter matrices
#' @param prev_naa the NAA in the previous timestep
#' @param recruitment total recruiment to occurin this year
#' @param option additional model options
#'
#' @return list of derived quantities included landed catch-at-age,
#' discarded catch-at-age, total catch-at-age, F-at-age, and numbers-at-age.
#'
#' @export
#'
project <- function(TAC, fleet.props, dem_params, prev_naa, recruitment, options=NA){

    model_params <- get_model_dimensions(dem_params$sel)

    if(!("region_apportionment" %in% names(options))){
        options$region_apportionment <- rep(1/model_params$nregions, model_params$nregions)
    }

    if(model_params$nregions < 2){
        options$region_apportionment <- 1
    }

    # TODO: make this into a list and add function to update the whole list at once
    land_caa_tmp    = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    disc_caa_tmp    = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    caa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    faa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    zaa_tmp         = array(0, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))
    naa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))

    # Do recruitment here because there isn't regional recruitment
    # NOTE: need to think real carefully about how to do this in a modular way
    rec <- array(NA, dim=c(1, 1, model_params$nsexes, model_params$nregions))
    global.rec <- recruitment

    # If no user-definied regional apportionment for recruitment than regionally
    # apportion based on the regional prevalence of age-1 individuals in the previous
    # timestep.
    recruit_apportionment <- NA
    if(model_params$nregions > 1){
        recruit_apportionment <- options$region_apportionment
        if(!all(is.na(recruit_apportionment))){
            rec.props <- prev_naa[1,1,1,]/sum(prev_naa[1,1,1,])
            multi <- rmultinom(50, model_params$nregions, prob=rec.props)
            recruit_apportionment <- apply(multi, 1, mean)/sum(apply(multi, 1, mean))
        }
        rec[1,1,,] <- t(apply(global.rec * dem_params$sexrat[1,1,,], 1, \(x) x*recruit_apportionment))
    }else{
        recruit_apportionment <- 1
        rec[1,1,,] <- global.rec * dem_params$sexrat[1,1,,]
    }

    for(r in 1:model_params$nregions){
        tac <- TAC*options$region_apportionment[r]
        dp.r <- subset_dem_params(dem_params=dem_params, r=r, d=4, drop=FALSE)
        prev_naa <- subset_dem_params(prev_naa, r=r, d=4, drop=FALSE)
        catch_vars <- simulate_catch(TAC=tac, dem_params=dp.r, naa=prev_naa, fleet.props = fleet.props, options=options)

        land_caa_tmp[,,,r,] <- catch_vars$land_caa_tmp
        disc_caa_tmp[,,,r,] <- catch_vars$disc_caa_tmp
        caa_tmp[,,,r,] <- catch_vars$caa_tmp
        faa_tmp[,,,r,] <- catch_vars$faa_tmp

        tot_faa <- array(apply(faa_tmp, c(2, 3), sum), dim=c(1, model_params$nages, model_params$nsexes, 1))
        zaa_tmp[,,,r] <- tot_faa+dem_params$mort

        pop_vars <- simulate_population(prev.naa=prev_naa, faa=catch_vars$faa_tmp, recruitment=rec, dem_params=dp.r, options=options)
        naa_tmp[,,,r] <- pop_vars$naa
    }
    # state_vars <- simulate_movement(dem_params, state_vars)
    obs <- simulate_observations(naa_tmp, dem_params$waa, dem_params$surv_sel, faa_tmp, zaa_tmp, obs_pars = list(surv_ll_q=6.14, surv_tw_q=0.851, ll_idx=1))
    surv_preds <- obs$preds
    surv_obs <- obs$obs


    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, faa_tmp, naa_tmp, surv_preds, surv_obs))
}
