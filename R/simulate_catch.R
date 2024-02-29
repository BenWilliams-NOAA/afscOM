#' Simulate catches in one region during on year
#' 
#' Description
#'
#' @param TAC the total catch to be removed across all fleets
#' @param fleet.props the proportion of the TAC that should be taken by each fleet
#' @param dem_params demographic parameter matrices for population
#' @param naa the current numbers-at-age in the popultion
#' @param option list of model options
#'
#' @export simualte_catch
#'
#' @example
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
        tac <- TAC*fleet.props[f]
        
        F_f <- find_F(
            f_guess = 0.05, 
            naa     = naa,
            waa     = dem_params$waa,
            mort    = dem_params$mort,
            selex   = dem_params$sel[,,f],
            ret     = dem_params$ret[,,f],
            dmr     = dem_params$dmr[,,f],
            prov_catch = tac
        )

        ret_faa <- retained_F(F_f, dem_params$sel[,,f], dem_params$ret[,,f])
        disc_faa <- discard_F(dem_params$dmr[,,f], dem_params$sel[,,f], dem_params$ret[,,f])
        faa_tmp[,,,,f] <- ret_faa + disc_faa

        land_caa_tmp[,,,,f] <- catch_at_age(ret_faa, naa, dem_params$waa, dem_params$mort)
        disc_caa_tmp[,,,,f] <- catch_at_age(disc_faa, naa, dem_params$waa, dem_params$mort)
        caa_tmp[,,,,f] <- land_caa_tmp[,,,,f] + disc_caa_tmp[,,,,f]
    }

    return(listN(land_caa_tmp, disc_caa_tmp, caa_tmp, faa_tmp))

}
