#' Generate a filled parameter matrix
#'
#' Fills an empty matrix of the dimensions (nyears, nages, nsexes, nregions, nfleets)
#' with a vector or matrix of values across the specified dimensions. This allows for
#' quickly filling parameter matrices when parameters only vary across a subset of the
#' required dimensions.
#'
#' @param vals a single value, or a vector or matrix of values, to fill the output matrix with. If a vector or matrix, all dimensions must be named, and names must match with those specified with the `dimension.names` argument
#' @param settings names list containing `model.params` that
#' @param by a character vector specifying which dimensions the input values correspond to
#' @param include_fleet_dim whether to expand the output matrix to include a 5th dimension, indicating fleet structure
#'
#' @return a 4d or 5d array of dimensions (nyears, nages, nsexes, nregions, nfleets) filled across the specified dimensions by the specfied values.
#'
#' @export
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
