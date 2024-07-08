#' Project Population Forward Multiple Years
#' 
#' A wrapper function around the base `project` function that handles
#' projecting forward multiple years given a removals timeseries and a 
#' recruitment timeseries.
#'
#' @param init_naa numbers-at-age matrix in starting year ([1, nages, nsexes, nregions])
#' @param removals_timeseries vector of removals (either catch of F) of length nyears
#' @param recruitment vector of recruitmene of length nyears
#' @param dem_params list of demographic parameter matrices
#' @param nyears number of projection yeas
#' @param model_option list of additional model options
#'
#' @export project_multi
#'
#' @example
#'
project_multi <- function(init_naa, removals_timeseries, recruitment, dem_params, nyears, model_options){

    model_dimensions <- get_model_dimensions(dem_params$sel)
    nyears <- model_dimensions$nyears
    nages <- model_dimensions$nages
    nsexes <- model_dimensions$nsexes
    nregions <- model_dimensions$nregions
    # nsurveys <- get_model_dimensions(dem_params$surv_sel)$nsurveys
    nsurveys <- ifelse(model_options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
    naa[1,,,] = init_naa

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))
    recruits    = array(NA, dim=c(nyears, 1, 1, nregions))

    survey_preds <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    survey_obs <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    # full_recruitment <- array(NA, dim=c(nyears, 1, 1, nregions))
    if(!(is.null(model_options$recruit_apportionment)) & !is.function(model_options$recruit_apportionment)){
        rec_props <- model_options$recruit_apportionment
    }else{
        rec_props <- array(1/nregions, dim=c(nyears+1, nregions))
    }

    # if recruitment is entered as a vector of global recruitment and
    # regional recruitment apportionment
    if(is.vector(recruitment) | is.array(recruitment) & dim(recruitment)[2] == 1){
        full_recruitment <- sweep(rec_props, 1, recruitment, FUN="*")
    }else if(all(is.array(recruitment) & dim(recruitment) > 1)){
        full_recruitment <- recruitment
    }


    for(y in 1:nyears){

        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- subset_matrix(removals_timeseries, y, d=1, drop=FALSE)
        # fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))
        # region_props <- model_options$region_apportionment[y,,drop=FALSE]
        # fleet_props <- model_options$fleet_apportionment[y,,drop=FALSE]
        region_props <- subset_matrix(model_options$region_apportionment, y, d=1, drop=FALSE)
        fleet_props <- subset_matrix(model_options$fleet_apportionment, y, d=1, drop=FALSE)

        if(!is.function(model_options$recruit_apportionment)){
            rec <- subset_matrix(full_recruitment, y+1, d=1, drop=FALSE)
        }else {
            projected_rec_props <- do.call(
                model_options$recruit_apportionment, 
                c(list(naa=naa[y,,,,drop=FALSE], dem_params=dp.y), model_options$recruitment_pars)
            )
            rec_props <- array(projected_rec_props, dim=c(1, nregions))
            rec <- array(recruitment[y+1]*projected_rec_props, dim=c(1, nregions))
        }
        
        if(!is.null(model_options$recruit_apportionment_random) & model_options$recruit_apportionment_random){
            rand_rec_props <- rmultinom(1, size=30, prob = rec_props[y+1,])
            rand_rec_props <- rand_rec_props/sum(rand_rec_props)
            rec <- array(recruitment[y+1]*rand_rec_props, dim=c(1, nregions))
        }


        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=naa[y,,,, drop = FALSE],
            recruitment=rec,
            region_props = region_props,
            fleet_props = fleet_props,
            # rec_props = rec_props,
            options=model_options
        )

        # update state
        land_caa[y,,,,] <- out_vars$land_caa_tmp
        disc_caa[y,,,,] <- out_vars$disc_caa_tmp
        caa[y,,,,] <- out_vars$caa_tmp
        faa[y,,,,] <- out_vars$faa_tmp
        naa[y+1,,,] <- out_vars$naa_tmp

        f[y,,,,] <- out_vars$F_f_tmp
        recruits[y,,,] <- rec

        if(model_options$simulate_observations){
            survey_preds$rpns[y,,,,] <- out_vars$survey_preds$rpns
            survey_preds$rpws[y,,,,] <- out_vars$survey_preds$rpws
            survey_preds$acs[y,,,,]  <- out_vars$survey_preds$acs

            survey_obs$rpns[y,,,,] <- out_vars$survey_obs$rpns
            survey_obs$rpws[y,,,,] <- out_vars$survey_obs$rpws
            survey_obs$acs[y,,,,]  <- out_vars$survey_obs$acs
        }

    }

    return(listN(land_caa, disc_caa, caa, faa, f, naa, recruits, survey_preds, survey_obs))

}
