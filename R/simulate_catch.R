#' Simulate catches in one region during on year
#'
#' Description
#'
#' @param TAC the total catch to be removed across all fleets
#' @param fleet.props the proportion of the TAC that should be taken by each fleet
#' @param dem_params demographic parameter matrices for population
#' @param naa the current numbers-at-age in the population
#' @param option list of model options
#'
#' @export
#'
simulate_catch <- function(TAC, fleet.props, dem_params, naa, options){

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
        if(model_options$fleet_is_survey[f]){
            next;
        }

        tac <- TAC*fleet.props[f]

        dp.f <- rlang::duplicate(dem_params)
        dp.f$sel <- subset_matrix(dp.f$sel, f, d=5, drop=TRUE)
        dp.f$ret <- subset_matrix(dp.f$ret, f, d=5, drop=TRUE)
        dp.f$dmr <- subset_matrix(dp.f$dmr, f, d=5, drop=TRUE)

        F_f <- findF_bisection(
            f_guess = 0.05,
            naa     = naa,
            waa     = dp.f$waa,
            mort    = dp.f$mort,
            selex   = dp.f$sel,
            ret     = dp.f$ret,
            dmr     = dp.f$dmr,
            prov_catch = tac
        )

        # ret_faa <- retained_F(F_f, dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        # disc_faa <- discard_F(dem_params$dmr[,,,,f], dem_params$sel[,,,,f], dem_params$ret[,,,,f])
        ret_faa <- retained_F(F_f, dp.f$sel, dp.f$ret)
        disc_faa <- discard_F(dp.f$dmr, dp.f$sel, dp.f$ret)
        faa_tmp[,,,,f] <- ret_faa + disc_faa

        land_caa_tmp[,,,,f] <- catch_at_age(ret_faa, naa, dp.f$waa, dp.f$mort)
        disc_caa_tmp[,,,,f] <- catch_at_age(disc_faa, naa, dp.f$waa, dp.f$mort)
        caa_tmp[,,,,f] <- land_caa_tmp[,,,,f] + disc_caa_tmp[,,,,f]
    }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, faa_tmp))

}
