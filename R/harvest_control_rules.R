constant_catch_hcr <- function(c, ...){
    return(c)
}

constant_f_hcr <- function(f, ...){
    return(f)
}

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


apply_harvest_control_rule <- function(model_dimensions, hcr_func, hcr_pars){
    expected_dim <- c(1, model_dimensions$nfleets, model_dimensions$nregions)
    hcr_out <- do.call(hcr_func, hcr_pars)
    if(!is.null(dim(hcr_out)) && all(dim(hcr_out) == expected_dim)){
        return(hcr_out)
    }else{
        # If have to coerce into proper output format, breakup catch/F so thats its
        # equivalent across regions and fleets. If more nuanced apportionments is
        # required, users mus handle that themselves.
        hcr_out_array <- array(hcr_out/model_dimensions$nregions/model_dimensions$nfleets, dim=expected_dim)
        return(hcr_out_array)
    }
}
