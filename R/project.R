#' Project Population Forward Multiple Years
#'
#' A wrapper function around the base `project` function that handles
#' projecting forward multiple years given a removals timeseries and a
#' recruitment timeseries.
#'
#' @param init_naa numbers-at-age matrix in starting year ([1, nages, nsexes, nregions])
#' @param removals_timeseries vector of removals (either catch of F) of length nyears
#' @param recruitment vector of recruitment of length nyears
#' @param dem_params list of demographic parameter matrices
#' @param nyears number of projection yeas
#' @param model_option list of additional model options
#'
#' @export project
#'
project <- function(init_naa, removals_timeseries, recruitment, dem_params, nyears, model_options){

    model_dimensions <- get_model_dimensions(dem_params$sel)
    nyears <- model_dimensions$nyears
    nages <- model_dimensions$nages
    nsexes <- model_dimensions$nsexes
    nregions <- model_dimensions$nregions
    nfleets <- model_dimensions$nfleets
    nsurveys <- ifelse(model_options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)

    outputs <- generate_output_matrices(nyears, nages, nsexes, nregions, nfleets, nsurveys)
    outputs$naa[1,,,] <- init_naa

    r <- apportion_recruitment(
        rec_timeseries = recruitment, 
        apportionment = model_options$recruit_apportionment,
        nyears = nyears,
        nregions = nregions
    )

    c <- apportion_catch(
        catch_timeseries = removals_timeseries,
        apportionment = model_options$fleet_apportionment,
        nyears = nyears,
        nfleets = nfleets,
        nregions = nregions
    )$full_catch

    set.seed(model_options$seed)
    for(y in 1:nyears){

        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        
        # Apply harvest control rule if insufficient catch/F timeseries is
        # provided as input...
        #
        # TODO: Handle case where no HCR supplied and insufficient catch/F timeseries
        if(y > dim(c)[1] && !is.null(model_options$hcr)){
            hcr_out <- apply_harvest_control_rule(model_dimensions, hcr_func=model_options$hcr$hcr_func, hcr_pars=c(model_options$hcr$hcr_pars, list(naa=outputs$naa[y,,,,drop=FALSE], dem_params=dp.y)))
            new_c <- array(NA, dim=c(dim(c)[1]+1, nfleets, nregions))
            new_c[1:dim(c)[1],,] <- c
            new_c[y,,] <- hcr_out
            c <- new_c
        }
        removals_input <- subset_matrix(c, y, d=1, drop=FALSE)
        
        

        if(is.function(recruitment)){
            r_y <- do.call(recruitment, c(list(naa=outputs$naa[y,,,,drop=FALSE], dem_params=dp.y), model_options$recruitment_pars))
        }else{
            rs <- array(0, dim=c(nyears+1, nregions))
            rs[1:nyears, 1:nregions] <- recruitment
            # rs <- array(recruitment, dim=c(nyears, nregions))
            # rs[y+1,] <- rep(0, nregions)
            r_y <- subset_matrix(rs, y+1, d=1, drop=FALSE)
        }
        r_y <- as.vector(r_y)

        r <- apportion_recruitment_single(
            recruits = as.vector(r_y),
            apportionment = subset_matrix(model_options$recruit_apportionment, y+1, d=1, drop=FALSE),
            nregions = nregions
        )

        rec <- get_annual_recruitment(
            recruitment = r$full_recruitment,
            apportionment = r$rec_props,
            apportion_random = model_options$recruit_apportionment_random,
            apportionment_pars = model_options$recruit_apportionment_pars,
            nregions = nregions,
            list(naa=outputs$naa[y,,,,drop=FALSE], dem_params=dp.y)
        )

        out_vars <- project_single(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=outputs$naa[y,,,, drop = FALSE],
            recruitment=rec,
            options=model_options
        )

        # update state
        outputs <- update_output_matrices(outputs, y, out_vars, update_obs = model_options$simulate_observations)
        outputs$recruits[y+1,,,] <- rec

    }

    return(outputs)
    #return(listN(land_caa, disc_caa, caa, faa, f, naa, recruits, survey_preds, survey_obs))

}
