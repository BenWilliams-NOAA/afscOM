#' Constant Catch Harvest Control Rule
#' 
#' Return a constant level of harvest (c) in
#' catch units.
#' 
#' @param c constant catch level
#' 
#' @return c
#' 
#' @export constant_catch_hcr 
#'
constant_catch_hcr <- function(c, ...){
    return(c)
}

#' Constant Fishing Mortality Harvest Control Rule
#' 
#' Return a constant level of fishing mortality (F) in
#' catch units.
#' 
#' @param f constant catch level
#' 
#' @return f
#' 
#' @export constant_f_hcr
#' 
constant_f_hcr <- function(f, ...){
    return(f)
}

#' Threshold Fishing Mortality Harvest Control Rule
#' 
#' Return an allowable fishing mortality rate based on
#' a threshold (hockey-stick) harvest control rule, using 
#' spawning stock biomass (SSB) as the input.
#' 
#' @param fmax maximum F to apply
#' @param fmin minimum F to apply
#' @param urp SSB level above which fmax is applied
#' @param lrp SSB level below which fmin is applied
#' @param naa numbers-at-age vector (used to internally compute SSB)
#' @param dem_params demographic parameters (used to internally compute SSB)
#' 
#' @return F as defined by SSB and the provided reference points
#' 
#' @export threshold_f_hcr 
#' 
threshold_f_hcr <- function(fmax, fmin, urp, lrp, naa, dem_params){
    ssb <- as.numeric(compute_ssb(naa, dem_params))
    x <- ssb
    if(x < lrp){
        return(fmin)
    }else if(x >= urp){
        return(fmax)
    }else{
        return(fmax*(x-lrp)/(urp-lrp)+fmin)
    }
}


#' Apply a harvest control rule
#'
#' @param model_dimensions a named list indicating the dimensions of the model
#' @param hcr_func the type of control rule
#' @param hcr_pars the associated parameters
#'
#' @export apply_harvest_control_rule
#'
apply_harvest_control_rule <- function(model_dimensions, hcr_func, hcr_pars){
    expected_dim <- c(1, model_dimensions$nfleets, model_dimensions$nregions)
    hcr_out <- do.call(hcr_func, hcr_pars)
    if(!is.null(dim(hcr_out)) && all(dim(hcr_out) == expected_dim)){
        return(hcr_out)
    }else{
        # If have to coerce into proper output format, breakup catch/F so that it's
        # equivalent across regions and fleets. If a more nuanced apportionment is
        # required, users must handle that themselves.
        hcr_out_array <- array(hcr_out/model_dimensions$nregions/model_dimensions$nfleets, dim=expected_dim)
        return(hcr_out_array)
    }
}
