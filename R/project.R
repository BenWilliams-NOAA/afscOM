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
project <- function(removals, dem_params, prev_naa, recruitment, region_props, fleet_props, rec_props, options=NA){

    model_params <- get_model_dimensions(dem_params$sel)

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

    # Do recruitment here because there isn't regional recruitment
    # NOTE: need to think real carefully about how to do this in a modular way
    rec <- array(NA, dim=c(1, 1, model_params$nsexes, model_params$nregions))
    global.rec <- recruitment

    # If no user-definied regional apportionment for recruitment than regionally
    # apportion based on the regional prevalence of age-1 individuals in the previous
    # timestep.
    recruit_apportionment <- NA
    if(model_params$nregions > 1){
        if(!(is.null(rec_props) | any(is.nan(rec_props)) | any(is.na(rec_props)))){
            recruit_apportionment <- rec_props
        }else{
            recruit_apportionment <- 1/model_params$nregions
        }
        # if(!all(is.na(recruit_apportionment))){
        #     rec.props <- prev_naa[1,1,1,]/sum(prev_naa[1,1,1,])
        #     multi <- rmultinom(50, model_params$nregions, prob=rec.props)
        #     recruit_apportionment <- apply(multi, 1, mean)/sum(apply(multi, 1, mean))
        # }
        rec[1,1,,] <- t(apply(global.rec * dem_params$sexrat[1,1,,], 1, \(x) x*recruit_apportionment))
    }else{
        recruit_apportionment <- 1
        rec[1,1,,] <- global.rec * dem_params$sexrat[1,1,,]
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
        catch_vars <- simulate_catch(
            removals=remove,
            dem_params=dp.r,
            naa=prev_naa.r,
            fleet_props = fleet_props,
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
    }

     # Handle movement matrix
    if(model_params$nregions > 1 & "movement" %in% names(dem_params)){
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
            FUN.VALUE = array(0, dim=c(model_params$nsexes, model_params$nregions))
        )
        moved_naa <- array(aperm(v, perm=c(3, 1, 2)), dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))
        naa_tmp <- moved_naa
    }

    # state_vars <- simulate_movement(dem_params, state_vars)
    surv_preds <- list()
    surv_obs <- list()

    if(!("simulate_observations" %in% names(options)) || options$simulate_observations){
        obs <- simulate_observations(
            naa = prev_naa,
            waa = dem_params$waa,
            selex = dem_params$surv_sel,
            faa = faa_tmp,
            zaa = zaa_tmp,
            obs_pars = options$obs_pars
        )
        surv_preds <- obs$preds
        surv_obs <- obs$obs
    }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, F_f_tmp, faa_tmp, naa_tmp, surv_preds, surv_obs))
}
