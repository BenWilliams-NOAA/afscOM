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
    # nsurveys <- get_model_dimensions(dem_params$surv_sel)$nsurveys
    nsurveys <- ifelse(model_options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
    naa[1,,,] = init_naa

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))
    recruits    = array(NA, dim=c(nyears+1, 1, 1, nregions))

    survey_preds <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    survey_obs <- list(
        catch = array(NA, dim=c(nyears, 1, 1, nregions, nfleets)), 
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    # full_recruitment <- array(NA, dim=c(nyears, 1, 1, nregions))
    # r <- apportion_recruitment(
    #     rec_timeseries = recruitment, 
    #     apportionment = model_options$recruit_apportionment,
    #     nyears = nyears,
    #     nregions = nregions
    # )

    if(model_options$removals_input == "catch"){
        c <- apportion_catch(
            catch_timeseries = removals_timeseries,
            apportionment = model_options$fleet_apportionment,
            nyears = nyears,
            nfleets = nfleets,
            nregions = nregions
        )$full_catch
    }else{
        c <- removals_timeseries
    }

    set.seed(model_options$seed)
    for(y in 1:nyears){

        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- subset_matrix(c, y, d=1, drop=FALSE)
        # removals_input <- subset_matrix(removals_timeseries, y, d=1, drop=FALSE)
        # region_props <- subset_matrix(model_options$region_apportionment, y, d=1, drop=FALSE)
        # fleet_props <- subset_matrix(model_options$fleet_apportionment, y, d=1, drop=FALSE)

        if(is.function(recruitment)){
            r_y <- do.call(recruitment, c(list(naa=naa[y,,,,drop=FALSE], dem_params=dp.y), model_options$recruitment_pars))
        }else{
            rs <- array(recruitment, dim=c(nyears+1, 1))
            r_y <- subset_matrix(rs, y+1, d=1, drop=FALSE)
        }

        # r_y <- ifelse(!is.null(model_options$recruitment_devs), as.vector(exp(log(r_y)+model_options$recruitment_devs[y+1])), as.vector(r_y))
        r_y <- as.vector(r_y)

        r <- apportion_recruitment_single(
            recruits = as.vector(r_y),
            apportionment = subset_matrix(model_options$recruit_apportionment, y+1, d=1, drop=FALSE),
            nregions = nregions
        )

        rec <- get_annual_recruitment(
            recruitment = r$full_recruitment,
            apportionment = r$rec_props,
            apportion_random = model_options$random_apportion_recruits,
            apportionment_pars = model_options$recruit_apportionment_pars,
            nregions = nregions,
            list(naa=naa[y,,,,drop=FALSE], dem_params=dp.y)
        )
        
        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=naa[y,,,, drop = FALSE],
            recruitment=rec,
            # region_props = region_props,
            # fleet_props = fleet_props,
            options=model_options
        )

        # update state
        land_caa[y,,,,] <- out_vars$land_caa_tmp
        disc_caa[y,,,,] <- out_vars$disc_caa_tmp
        caa[y,,,,] <- out_vars$caa_tmp
        faa[y,,,,] <- out_vars$faa_tmp
        naa[y+1,,,] <- out_vars$naa_tmp

        f[y,,,,] <- out_vars$F_f_tmp
        recruits[y+1,,,] <- rec

        if(model_options$simulate_observations){
            survey_preds$rpns[y,,,,] <- out_vars$survey_preds$rpns
            survey_preds$rpws[y,,,,] <- out_vars$survey_preds$rpws
            survey_preds$acs[y,,,,]  <- out_vars$survey_preds$acs

            survey_obs$catch[y,,,,] <- out_vars$survey_obs$catch
            survey_obs$rpns[y,,,,] <- out_vars$survey_obs$rpns
            survey_obs$rpws[y,,,,] <- out_vars$survey_obs$rpws
            survey_obs$acs[y,,,,]  <- out_vars$survey_obs$acs
        }

    }

    return(listN(land_caa, disc_caa, caa, faa, f, naa, recruits, survey_preds, survey_obs))

}
