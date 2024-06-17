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
    nfleets <- model_dimensions$nfleets

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
    naa[1,,,] = init_naa

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))

    survey_preds <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_ac = array(NA, dim=c(nyears, nages, nsexes, nregions)),
        fxfish_caa = array(NA, dim=c(nyears, nages, nsexes, nregions))
    )

    survey_obs <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions)),
        fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions))
    )


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
        rec_props <- subset_matrix(model_options$recruit_apportionment, y+1, d=1, drop=FALSE)

        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=naa[y,,,, drop = FALSE],
            recruitment=recruitment[y+1],
            region_props = region_props,
            fleet_props = fleet_props,
            rec_props = rec_props,
            options=model_options
        )

        # update state
        land_caa[y,,,,] <- out_vars$land_caa_tmp
        disc_caa[y,,,,] <- out_vars$disc_caa_tmp
        caa[y,,,,] <- out_vars$caa_tmp
        faa[y,,,,] <- out_vars$faa_tmp
        naa[y+1,,,] <- out_vars$naa_tmp

        f[y,,,,] <- out_vars$F_f_tmp

        if(model_options$simulate_observations){
            survey_preds$ll_rpn[y,,,] <- out_vars$surv_preds$ll_rpn
            survey_preds$ll_rpw[y,,,] <- out_vars$surv_preds$ll_rpw
            survey_preds$tw_rpw[y,,,] <- out_vars$surv_preds$tw_rpw
            survey_preds$ll_ac[y,,,] <- out_vars$surv_preds$ll_ac
            survey_preds$fxfish_caa[y,,,] <- out_vars$surv_preds$fxfish_caa

            survey_obs$ll_rpn[y,,,] <- out_vars$surv_obs$ll_rpn
            survey_obs$ll_rpw[y,,,] <- out_vars$surv_obs$ll_rpw
            survey_obs$tw_rpw[y,,,] <- out_vars$surv_obs$tw_rpw
            survey_obs$ll_acs[y,,,] <- out_vars$surv_obs$ll_ac_obs
            survey_obs$fxfish_acs[y,,,] <- out_vars$surv_obs$fxfish_caa_obs
        }

    }

    return(listN(land_caa, disc_caa, caa, faa, f, naa, survey_preds, survey_obs))

}
