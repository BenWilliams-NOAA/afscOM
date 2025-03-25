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

apportion_recruitment_single <- function(recruits, apportionment, nregions){
    
    rec_props <- NULL
    if(is.null(apportionment)){
        rec_props <- array(1/nregions, dim=c(1, nregions))
    }else if(is.vector(apportionment)){
        rec_props <- array(apportionment, dim=c(1, nregions))
    }else{
        rec_props <- apportionment
    }

    # if recruitment is entered as a vector of global recruitment and
    # regional recruitment apportionment
    full_recruitment <- recruits
    if(!is.function(rec_props) && !(length(full_recruitment) > 1)){
        full_recruitment <- recruits*rec_props
    }
    full_recruitment <- array(full_recruitment, dim=c(1, nregions))

    return(listN(rec_props, full_recruitment))
}

get_annual_recruitment <- function(recruitment, apportionment, apportion_random, apportionment_pars, nregions, ...){
    
    rec_props <- array(NA, dim=c(1, nregions))
    if(!is.function(apportionment)){
        rec <- recruitment
        rec_props <- rec/sum(rec)
    }else {
        projected_rec_props <- do.call(
            apportionment, 
            c(list(), apportionment_pars)
        )
        rec_props <- array(projected_rec_props, dim=c(1, nregions))
        rec <- array(recruitment*projected_rec_props, dim=c(1, nregions))
    }
    
    if(apportion_random){
        rand_rec_props <- rmultinom(1, size=30, prob = rec_props)
        rand_rec_props <- t(rand_rec_props/sum(rand_rec_props))
        rec <- array(recruitment*rand_rec_props, dim=c(1, nregions))
    }

    return(rec)

}

apportion_catch <- function(catch_timeseries, apportionment, nyears, nfleets, nregions){
    
    # Derive the recruitment proportions from an input catch
    # matrix if provided
    nyears <- ifelse(is.array(catch_timeseries), dim(catch_timeseries)[1], length(catch_timeseries))

    if(length(dim(catch_timeseries)) == 3 && all(dim(catch_timeseries) == c(nyears, nfleets, nregions))){
        region_fleet_props <- vapply(
            1:nyears, 
            function(y){
                catch_y <- catch_timeseries[y,,]
                catch_y/sum(catch_y)
            },
            FUN.VALUE = array(0, dim=c(nfleets, nregions))
        )
        region_fleet_props <- aperm(region_fleet_props, c(3, 1, 2))
        full_catch <- catch_timeseries
        return(listN(region_fleet_props, full_catch))
    }

    region_fleet_props <- NULL
    if(is.null(apportionment)){
        region_fleet_props <- array(1/nregions, dim=c(nyears, nfleets, nregions))
    }else if(is.vector(apportionment)){
        region_fleet_props <- array(apportionment, dim=c(nyears, nfleets, nregions))
    }else{
        region_fleet_props <- apportionment[1:nyears,,,drop=FALSE]
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
