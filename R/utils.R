#' harvest rate to instantaneous
#'
#' Converts between harvest rate - the proportion
#' of the population harvest; mu and instantaneous
#' fishing mortality rate F_
#'
#' @param mu, a harvest rate
#'
#' @return F, the corresponding instantaneous fishing mortality rate
#'
#' @export mu_to_F
#'
mu_to_F <- function(mu){
    return(-log(1-mu))
}

#' instantaneous rate to harvest rate
#'
#' Converts between instantaneous fishing mortality
#' rate F and harvest rate -the proportion of the
#' population harvest; mu_
#'
#' @param fy, an instantaneous fishing mortality rate
#'
#' @return mu, the corresponding harvest rate
#'
#' @export F_to_mu
#'
F_to_mu <- function(fy){
    return(1-exp(-fy))
}

#' Set variable names as names of list
#'
#' Automatically use variable names as list object names when
#' a sequence of variables is provided to a list constructor_
#'
#' @param ... R variables containung values to be put in a list
#'
#' @export listN
listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

#' Subset demographic parameter matrices by first dimensions
#'
#' Access an arbitrary row across all demographic parameter
#' matrices of indeterminate dimension
#'
#' @param dem_params a list of arrays with more than 2 dimensions
#' @param r the index along the first dimension to access
#'
#' @return a list with the same elements as dem_params but containing
#' a single row from each list element
#'
#' @export subset_dem_params
#'
#'
subset_dem_params <- function(dem_params, r, d=1, drop=TRUE){
    tmp <- rlang::duplicate(dem_params)
    ps <- names(tmp)
    for(n in ps){
        if(n == "movement") next;
        tmp[[n]] <- subset_matrix(dem_params[[n]], r, d, drop)
    }
    return(tmp)
}

#' Subset matrix of undetermined dimensions
#'
#' Return a subset or a multidimensional matrix or array of undetermined
#' size along a specific dimension. Specifically designed to return a
#' single index along the specified dimension (e.g. one year or one region).
#'
#' @param mat the input matrix to subset
#' @param r the index along the dimensions of interest to subset by
#' @param d the dimension to subset by
#' @param drop whether to drop the subsetted dimension
#'
#' @return a subsetted matrix or array
#'
#' @export subset_matrix
#'
#'
subset_matrix <- function(mat, r, d=1, drop=TRUE){
    tmp <- rlang::duplicate(mat)
    ndims <- length(dim(mat))
    idxs <- c(as.list(rep(TRUE, d-1)), list(r), as.list(rep(TRUE, ndims-d)))
    t <- do.call('[', c(list(mat), idxs, drop=FALSE))
    if(drop){
        tmp <- abind::adrop(t, drop=d)
    }else{
        tmp <- t
    }
    return(tmp)
}

# TODO: generalize this to allow extending across other dimensions

#' Extend demographic parameter matrices to more years
#'
#' Description
#'
#' @param dem_params a list of demographic parameter matrices
#' @param dimension the dimension along which to extend (should always be 1)
#' @param e the new number of years along the first dimension
#' @param new.dimnames new set of names for the first dimension
#'
#' @return a new list of demographic parameter matrices spanning `e` years
#' with the values from the last year in the original matrix continued for
#' all future years
#'
#' @export extend_years
#'
#'
extend_years <- function(dem_params, dimension, e, new.dimnames=NA){
    tmp <- rlang::duplicate(dem_params)
    ps <- names(tmp)
    for(n in ps){
        ndims <- length(dim(dem_params[[n]]))
        new.dims <- dim(dem_params[[n]])
        new.dims[dimension] <- e
        new.dimnames <- if(all(is.na(new.dimnames))) 1:e else new.dimnames
        t <- array(NA, dim=new.dims, dimnames = c("time"=list(new.dimnames), dimnames(dem_params[[n]])[2:length(dimnames(dem_params[[n]]))]))
        last <- subset_matrix(dem_params[[n]], r=nrow(dem_params[[n]]), d=1)
        afill(t) <- dem_params[[n]]

        afill_dimensions <- as.vector(c(TRUE, rep(FALSE, ndims-1)), mode="list")
        for(i in which(afill_dimensions == FALSE)){
            afill_dimensions[[i]] <- rlang::missing_arg()
        }
        afill_params <- afill_dimensions
        afill_params$x <- t
        afill_params$value <- last

        t <- do.call("afill<-", afill_params)
        tmp[[n]] <- t
    }
    return(tmp)
}


#' Set Default Values for Model Options
#' 
#' Set up a fully formed model_options list object with all
#' required elements set to sensible default values.
#'
#' @param model_dimensions model dimensions list object
#'
#' @export setup_model_options
#'
#' @example \dontrun{
#'      dimensions = list(nyears=60, nages=50, nsexes=2, nregions=1, nfleets=1)
#'      model_options = setup_model_options(dimensions)
#' }
#'
setup_model_options <- function(model_dimensions){

    return(
        list(
            removals_input = "catch",
            simulate_observations = TRUE,
            region_apportionment = list(rep(1/model_dimensions$nregions, model_dimensions$nregions)),
            fleet_apportionment = list(rep(1/model_dimensions$nfleets, model_dimensions$nfleets)),
            obs_pars = list(
                surv_ll = list(
                    q = 1,
                    rpn_cv = 0.10,
                    rpw_cv = 0.10,
                    ac_samps = 100,
                    as_integers = TRUE
                ),
                surv_tw = list(
                    q = 1,
                    rpw_cv = 0.10,
                    ac_samps = 100,
                    as_integers = TRUE
                ),
                fish_fx = list(
                    ac_samps = 100,
                    as_integers = TRUE
                ),
                fish_tw = list(
                    ac_samps = 100,
                    as_integers = TRUE
                )
            )
        )
    )
}