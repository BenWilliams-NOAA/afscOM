#' Project population forward by 1 year
#'
#' Given a TAC, the previous population age structure, demographic information
#' (e.g. natural mortality, selectivity, retention, etc.), and recruitment,
#' determine the numbers-at-age present in the population for the following year.
#'
#' @param removals population removals (either in catch or F units)
#' @param fleet.props the proportion of the TAC allocated to each fleet
#' @param dem_params list of demographic parameter matrices
#' @param prev_naa the NAA in the previous timestep
#' @param recruitment total recruiment to occurin this year
#' @param option additional model options
#'
#' @return list of derived quantities included landed catch-at-age,
#' discarded catch-at-age, total catch-at-age, F-at-age, and numbers-at-age.
#'
#' @export project
#'
project <- function(removals, dem_params, prev_naa, recruitment, region_props, fleet_props, options=NA){

    model_params <- get_model_dimensions(dem_params$sel)
    model_params$nsurveys <- ifelse(options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)
    # if(!("region_apportionment" %in% names(options))){
    #     options$region_apportionment <- rep(1/model_params$nregions, model_params$nregions)
    # }

    # if(model_params$nregions < 2){
    #     options$region_apportionment <- 1
    # }

    # TODO: make this into a list and add function to update the whole list at once
    land_caa_tmp    = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    disc_caa_tmp    = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    caa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    faa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
    F_f_tmp         = array(NA, dim=c(1, 1, 1,  model_params$nregions, model_params$nfleets))
    zaa_tmp         = array(0,  dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))
    naa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))

    survey_preds <- list(
        rpns = array(NA, dim=c(1, 1, 1, model_params$nregions, model_params$nsurveys)),
        rpws = array(NA, dim=c(1, 1, 1, model_params$nregions, model_params$nsurveys)),
        acs  = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nsurveys+model_params$nfleets))
    )

    survey_obs <- list(
        rpns = array(NA, dim=c(1, 1, 1, model_params$nregions, model_params$nsurveys)),
        rpws = array(NA, dim=c(1, 1, 1, model_params$nregions, model_params$nsurveys)),
        acs  = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nsurveys+model_params$nfleets))
    )

    # Do recruitment here because there isn't regional recruitment
    # NOTE: need to think real carefully about how to do this in a modular way
    rec <- array(NA, dim=c(1, 1, model_params$nsexes, model_params$nregions))
    global.rec <- as.vector(recruitment)
    if(model_params$nregions > 1){
        rec[1,1,,] <- sweep(dem_params$sexrat[1,1,,], 2, global.rec, FUN="*")
    }else{
        rec[1,1,,] <- global.rec*dem_params$sexrat[1,1,,]
    }

    for(r in 1:model_params$nregions){

        if(options$removals_input == "catch"){
            # Apportion catch-based removals based on provided
            # regional apportionment scheme.
            remove <- removals*region_props[1,r]
        }else{
            # Removals were input as F, subset to correct dimensions
            remove <- subset_matrix(removals, r=r, d=4, drop=FALSE)
        }

        dp.r <- subset_dem_params(dem_params=dem_params, r=r, d=4, drop=FALSE)
        prev_naa.r <- subset_matrix(prev_naa, r=r, d=4, drop=FALSE)

        if(length(dim(fleet_props)) > 2){
            fleet_props.r <- subset_matrix(fleet_props, r, d=3, drop=TRUE)
        }else{
            fleet_props.r <- fleet_props
        }

        catch_vars <- simulate_catch(
            removals=remove,
            dem_params=dp.r,
            naa=prev_naa.r,
            fleet_props = fleet_props.r,
            options=options
        )

        land_caa_tmp[,,,r,] <- catch_vars$land_caa_tmp
        disc_caa_tmp[,,,r,] <- catch_vars$disc_caa_tmp
        caa_tmp[,,,r,] <- catch_vars$caa_tmp
        faa_tmp[,,,r,] <- catch_vars$faa_tmp
        F_f_tmp[,,,r,] <- catch_vars$F_f


        tot_faa <- array(apply(catch_vars$faa_tmp, c(2, 3), sum), dim=c(1, model_params$nages, model_params$nsexes, 1))
        zaa_tmp[,,,r] <- tot_faa+dp.r$mort

        rec.r <- subset_matrix(rec, r=r, d=4, drop=FALSE)
        pop_vars <- simulate_population(
                        prev_naa=prev_naa.r, 
                        faa=catch_vars$faa_tmp, 
                        recruitment=rec.r, 
                        dem_params=dp.r, 
                        options=options
                    )
        naa_tmp[,,,r] <- pop_vars$naa

        if(!("simulate_observations" %in% names(options)) || options$simulate_observations){
            # concatenate fishery and survey selectivity matrices for convenience later
            big_selex <- abind::abind(dp.r$sel, dp.r$surv_sel, along=5) 
            names(dim(big_selex)) <- names(dim(dp.r$sel))
            obs <- simulate_observations(
                naa = prev_naa.r,
                waa = dp.r$waa,
                selex = big_selex,
                faa = faa_tmp[,,,r,,drop=FALSE],
                zaa = zaa_tmp[,,,r,drop=FALSE],
                obs_pars = options$obs_pars
            )
            survey_preds$rpns[,,,r,] <- obs$preds$rpn_preds[,,,as.logical(options$obs_pars$is_survey)]
            survey_preds$rpws[,,,r,] <- obs$preds$rpw_preds[,,,as.logical(options$obs_pars$is_survey)]
            survey_preds$acs[,,,r,]  <- obs$preds$ac_preds

            survey_obs$rpns[,,,r,] <- obs$obs$rpn_obs[,,,as.logical(options$obs_pars$is_survey)]
            survey_obs$rpws[,,,r,] <- obs$obs$rpw_obs[,,,as.logical(options$obs_pars$is_survey)]
            survey_obs$acs[,,,r,]  <- obs$obs$ac_obs
        }

    }

     # Handle movement matrix
    if(model_params$nregions > 1 & "movement" %in% names(dem_params)){
        if(!options$do_recruits_move){
            dem_params$movement[,,1,] <- diag(nregions)
        }
       v <- vapply(
            1:nages, 
            # Apply movement to the ages individualy
            function(a) {
                sapply(
                    1:model_params$nsexes,
                    # Apply movement to the sexes individually
                    function(s) naa_tmp[1,a,s,] %*% dem_params$movement[,,a,s]
                )
            } , 
            FUN.VALUE = array(0, dim=c(model_params$nregions, model_params$nsexes))
        )
        moved_naa <- array(aperm(v, perm=c(3, 2, 1)), dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))
        naa_tmp <- moved_naa
    }

    # state_vars <- simulate_movement(dem_params, state_vars)

    # if(!("simulate_observations" %in% names(options)) || options$simulate_observations){
    #     # concatenate fishery and survey selectivity matrices for convenience later
    #     big_selex <- abind::abind(dem_params$sel, dem_params$surv_sel, along=5) 
    #     names(dim(big_selex)) <- names(dim(dem_params$sel))
    #     obs <- simulate_observations(
    #         naa = prev_naa,
    #         waa = dem_params$waa,
    #         selex = big_selex,
    #         faa = faa_tmp,
    #         zaa = zaa_tmp,
    #         obs_pars = options$obs_pars
    #     )
    #     surv_preds <- obs$preds
    #     surv_obs <- obs$obs
    # }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, F_f_tmp, faa_tmp, naa_tmp, survey_preds, survey_obs))
}
