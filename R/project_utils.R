calculate_recruitment <- function(rec_timeseries, apportionment, nyears, nregions){
    
    rec_props <- NULL
    if(is.null(apportionment)){
        rec_props <- array(1/nregions, dim=c(nyears+1, nregions))
    }else if(is.vector(apportionment)){
        rec_props <- matrix(apportionment, nrow=nyears+1, ncol=nregions, byrow=TRUE)
    }else{
        rec_props <- apportionment
    }
    
    # if(!is.null(apportionment) & !is.function(apportionment)){
    #     rec_props <- apportionment
    # }

    # if recruitment is entered as a vector of global recruitment and
    # regional recruitment apportionment
    rec_timeseries <- array(rec_timeseries, dim=c(nyears+1, 1))
    full_recruitment <- rec_timeseries
    if(!is.function(rec_props)){
        full_recruitment <- sweep(rec_props, 1, rec_timeseries, FUN="*")
    }
    
    # if(is.vector(rec_timeseries) || is.array(rec_timeseries) & dim(rec_timeseries)[2] == 1){
        
    # }else if(all(is.array(rec_timeseries) & dim(rec_timeseries) > 1)){
    #     full_recruitment <- rec_timeseries
    # }

    return(listN(rec_props, full_recruitment))
}
