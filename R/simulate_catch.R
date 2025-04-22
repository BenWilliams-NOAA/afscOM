#' Simulate catches in one region during on year
#'
#' Description
#'
#' @param removals desired removals either in catch or F units
#' @param dem_params demographic parameter matrices for population
#' @param naa the current numbers-at-age in the population
#' @param options list of model options
#'
#' @export
#' @examples
#' /dontrun{
#' simulate_catch(removals, dem_params, naa, options)
#' }
simulate_catch <- function(removals, dem_params, naa, options){

    model_params <- get_model_dimensions(dem_params$sel)

    # if not provided, assume equal fleet split
    # if(!("fleet_apportionment" %in% names(options))){
    #     options$fleet_apportionment <- rep(1/model_params$nfleets, model_params$nfleets)
    # }

    land_caa_tmp    <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    disc_caa_tmp    <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    caa_tmp         <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    faa_tmp         <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1, model_params$nfleets))
    F_f_tmp         <- array(NA, dim=c(1, 1, 1, 1, model_params$nfleets))

    for(f in 1:model_params$nfleets){

        # Skip if the fleet in question is a survey fleet
        # if("fleet_is_survey" %in% model_options & model_options$fleet_is_survey[f]){
        #     next;
        # }

        dp.f <- rlang::duplicate(dem_params)
        dp.f$sel <- subset_matrix(dp.f$sel, f, d=5, drop=TRUE)
        dp.f$ret <- subset_matrix(dp.f$ret, f, d=5, drop=TRUE)
        dp.f$dmr <- subset_matrix(dp.f$dmr, f, d=5, drop=TRUE)

        if(options$removals_input == "catch"){
            # Apportion catch-based removals based on provided
            # fleet apportionment scheme.
            # remove <- removals[,f]
            remove <- as.numeric(subset_matrix(removals, f, d=2, drop=TRUE))

            # Solve for F that removes catch
            F_f <- find_F(
                f_guess = 0.05,
                naa     = naa,
                waa     = dp.f$waa,
                mort    = dp.f$mort,
                selex   = dp.f$sel,
                ret     = dp.f$ret,
                dmr     = dp.f$dmr,
                prov_catch = remove
            )
        }else{
            # Removals were input as F
          F_f <- removals[,f,]
        }

        # ret_faa <- retained_F(F_f, dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        # disc_faa <- discard_F(dem_params$dmr[,,,,f], dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        ret_faa <- retained_F(F_f, dp.f$sel, dp.f$ret)
        disc_faa <- discard_F(dp.f$dmr, dp.f$sel, dp.f$ret)
        faa_tmp[,,,,f] <- ret_faa + disc_faa

        land_caa_tmp[,,,,f] <- catch_at_age(ret_faa, naa, dp.f$waa, dp.f$mort)
        disc_caa_tmp[,,,,f] <- catch_at_age(disc_faa, naa, dp.f$waa, dp.f$mort)
        caa_tmp[,,,,f] <- land_caa_tmp[,,,,f] + disc_caa_tmp[,,,,f]
        F_f_tmp[,,,,f] <- F_f
    }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, faa_tmp, F_f_tmp))

}
