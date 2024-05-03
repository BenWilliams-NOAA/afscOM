#' Project an Age Structure Population Forward by 1 year
#' #'
#' Description
#'
#' @param prev_naa numbers-at-age vector for previous year
#' @param faa fishing mortality-at-age vector (dim [1, nages, nsexes, 1, nfleets])
#' @param recruitment projected recruitment for the folliwing year
#' @param dem_params demographic parameter matrices subsetted to 1 year and 1 region
#' @param options model options list
#'
#' @export simulate_population
#' @examples
#'\dontrun{
#'simulate_population(prev_naa, faa, recruitment, dem_params, options)}
#'
simulate_population <- function(prev_naa, faa, recruitment, dem_params, options){

    model_params <- get_model_dimensions(dem_params$sel)

    faa <- array(apply(faa, c(2, 3), sum), dim=c(1, model_params$nages, model_params$nsexes, 1)) # Sum across fleets to get total FAA
    Zaa <- faa + dem_params$mort

    naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))

    # Recruitment goes here
    naa[1,1,,1] <- recruitment

    # Simulate mortality here
    for(a in 2:(model_params$nages-1)){
        naa[1, a, ,1] <- prev_naa[1, a-1,,1]*exp(-Zaa[,a-1,,])
    }
    naa[1, model_params$nages, ,1] <- prev_naa[1, model_params$nages-1,,1]*exp(-Zaa[,model_params$nages-1,,])+prev_naa[1, model_params$nages,,1]*exp(-Zaa[,model_params$nages,,])

    return(listN(naa))

}
