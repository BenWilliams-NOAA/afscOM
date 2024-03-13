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

#' Generate a filled parameter matrix
#'
#' Fills an empty matrix of the dimensions (nyears, nages, nsexes, nregions, nfleets)
#' with a vector or matrix of values across the specified dimensions. This allows for
#' quicky filling parameter matrices when parameters only vary across a subset of the
#' required dimensions.
#'
#' @param vals a single value, or a vector or matrix of values, to fill the output matrix with. If a vector or matrix, all dimensions must be named, and names must match with those specified with the `dimension.names` argument
#' @param settings names list containing `model.params` that
#' @param by a character vector specifying which dimensions the input values correspond to
#' @param include_fleet_dim whether to expand the output matrix to include a 5th dimension, indicating fleet structure
#'
#' @return a 4d or 5d array of dimensions (nyears, nages, nsexes, nregions, nfleets) filled across the specified dimensions by the specfied values.
#'
#' @export generate_param_matrix
#'
#'

generate_param_matrix <- function(vals, dimension_names, by=NA, include_fleet_dim=FALSE){

    # if(any(is.na(model.params)) & !exists("model.params")){
    #     stop("`model.params` object does not exist and was not provided.")
    # }

    # if(any(is.na(dimension.names)) & !exists("dimension.names")){
    #     stop("`dimension.names` object does not exist and was not provided.")
    # }

    model_params <- lapply(dimension_names, length)
    names(model_params) <- c("nyears", "nages", "nsexes", "nregions", "nfleets")[1:length(model_params)]

    if("fleet" %in% by){
        include_fleet_dim <- TRUE
    }

    dim_order <- c("time", "age", "sex", "region", "fleet")

    ndims <- 4
    if(include_fleet_dim){
        ndims <- 5
    }

    if(length(vals)==1){
        arr <- array(vals, dim=model_params[1:ndims], dimnames=dimension_names[1:ndims])
        return(arr)
    }

    tmp <- array(NA, dim=model_params[1:ndims], dimnames=dimension_names[1:ndims])
    val_dims <- dim(vals)

    # vals is a c vector
    if(is.null(val_dims)){
        val_dims <- length(vals)
    }
    n_dims <- val_dims
    expected_dims <- unlist(lapply(dimension_names[by], length))
    if(!all(n_dims == expected_dims)){
        stop(paste0("Wrong dimensions (", n_dims ,") for filling by ", by, ". Expected dimensions are ", expected_dims, ".\n"))
    }
    vals <- array(vals, dim=val_dims, dimnames=dimension_names[by])
    provided_dimensions <- dim_order %in% by
    afill_dimensions <- as.vector(!provided_dimensions, mode="list")
    for(i in which(provided_dimensions == TRUE)){
        afill_dimensions[[i]] <- rlang::missing_arg()
    }
    afill_params <- afill_dimensions[1:ndims]
    afill_params$x <- tmp
    afill_params$value <- vals

    tmp <- do.call(abind::"afill<-", afill_params)
    return(tmp)
}

#' Obtain model dimensions from parameter matrix
#'
#' Automatically determines nyears, nages, nsexes, nregions, and nfleets
#' based on the dimensions of the input parameter matrix. This should be
#' run on the parameter matrix with the largest number of dimensions but
#' will return the correct dimensions for a smaller matrix.
#'
#' @param par_mat a parameter matrix generated via generate_param_matrix
#'
#' @return a named list indicating the dimensions of the model
#'
#' @export get_model_dimensions
#'
#'
get_model_dimensions <- function(par_mat){
    return(lapply(split(dim(par_mat), names(dim(par_mat))), unname))
}

set_model_params <- function(nyears, nages, nsexes=1, nregions=1, nfleets=1){
    return(
        list(
            nyears      = nyears,
            nages       = nages,
            nsexes      = nsexes,
            nregions    = nregions,
            nfleets     = nfleets
        )
    )
}