simulate_population <- function(prev.naa, faa, recruitment, dem_params, options){

    model_params <- get_model_dimensions(dem_params$sel)

    faa <- apply(faa, c(1,2), sum) # Sum across fleets to get total FAA
    Zaa <- faa + dem_params$mort

    naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))

    # Recruitment goes here
    naa[1,1,,1] <- recruitment

    # Simulate mortality here
    for(a in 2:(model_params$nages-1)){
        naa[1, a, ,1] <- prev.naa[1, a-1,,1]*exp(-Zaa[a-1,])
    }
    naa[1, model_params$nages, ,1] <- prev.naa[1, model_params$nages-1,,1]*exp(-Zaa[model_params$nages-1,])+prev.naa[1, model_params$nages,,1]*exp(-Zaa[model_params$nages,])

    return(listN(naa))

}