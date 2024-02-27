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
#' @example
#'
get_model_dimensions <- function(par_mat){
    return(lapply(split(dim(par_mat), names(dim(par_mat))), unname))
}