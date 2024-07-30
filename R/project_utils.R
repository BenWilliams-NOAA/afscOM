apportion_recruitment <- function(rec_timeseries, apportionment, nyears, nregions){
    
    rec_props <- NULL
    if(is.null(apportionment)){
        rec_props <- array(1/nregions, dim=c(nyears+1, nregions))
    }else if(is.vector(apportionment)){
        rec_props <- matrix(apportionment, nrow=nyears+1, ncol=nregions, byrow=TRUE)
    }else{
        rec_props <- apportionment
    }

    # if recruitment is entered as a vector of global recruitment and
    # regional recruitment apportionment
    rec_timeseries <- array(rec_timeseries, dim=c(nyears+1, 1))
    full_recruitment <- rec_timeseries
    if(!is.function(rec_props)){
        full_recruitment <- sweep(rec_props, 1, rec_timeseries, FUN="*")
    }

    return(listN(rec_props, full_recruitment))
}

get_annual_recruitment <- function(y, rec_timeseries, apportionment, apportion_random, apportionment_pars, nregions, ...){
    
    rec_props <- array(NA, dim=c(1, nregions))
    if(!is.function(apportionment)){
        rec <- subset_matrix(rec_timeseries, y+1, d=1, drop=FALSE)
        rec_props <- rec/sum(rec)
    }else {
        projected_rec_props <- do.call(
            apportionment, 
            c(list(), apportionment_pars)
        )
        rec_props <- array(projected_rec_props, dim=c(1, nregions))
        rec <- array(rec_timeseries[y+1,]*projected_rec_props, dim=c(1, nregions))
    }
    
    if(apportion_random){
        rand_rec_props <- rmultinom(1, size=30, prob = rec_props)
        rand_rec_props <- rand_rec_props/sum(rand_rec_props)
        rec <- array(rec_timeseries[y+1,]*rand_rec_props, dim=c(1, nregions))
    }

    return(rec)

}

apportion_catch <- function(catch_timeseries, apportionment, nyears, nfleets, nregions){
    
    region_fleet_props <- NULL
    if(is.null(apportionment)){
        region_fleet_props <- array(1/nregions, dim=c(nyears, nfleets, nregions))
    }else if(is.vector(apportionment)){
        region_fleet_props <- array(apportionment, dim=c(nyears, nfleets, nregions))
    }else{
        region_fleet_props <- apportionment
    }

    # if recruitment is entered as a vector of global recruitment and
    # regional recruitment apportionment
    catch_timeseries <- array(catch_timeseries, dim=c(nyears+1, 1))
    full_catch <- catch_timeseries
    if(!is.function(region_fleet_props)){
        full_catch <- sweep(region_fleet_props, 1, catch_timeseries, FUN="*")
    }

    return(listN(region_fleet_props, full_catch))
}
