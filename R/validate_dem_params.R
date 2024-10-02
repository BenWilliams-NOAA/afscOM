#' Validate Demographic Parameter Matrices
#'
#' Helper function for validating that all required demographic
#' parameter matrices have been defined and that they are all of
#' the correct dimension.
#'
#' Required demographic parameters are: mort, waa, sexrat, sel,
#' ret, dmr, and surv_sel
#'
#' Expected dimensions are: [nyears, nages, nsexes, nregions, nfleets]
#'
#' @param dem_params a named list containing the demographic rate matrices.
#' List elements must have the same names as stated above.
#' @param model_dims a named list containing expected model dimensions. List
#' element must be: nyears, nages, nsexes, nregions, and nfleets.
#'
#' @return unmodified dem_params object; if dem_params list is invalid, will
#' return a descriptive error message
#' @export
#'
#' @examples
#' \dontrun{
#'
#'      model_dims = list(nyears=64, nages=30, nsexes=2, nregions=1, nfleets=2)
#'      dem_params = list(
#'          mort = mortality_matrix,
#'          waa = weight_matrix,
#'          sexrat = sexratio_matrix,
#'          sel = selectivity_matrix,
#'          ret = retention_matrix
#'          dmr = disard_mortality_matrix,
#'          surv_sel = survey_selectivity_matrix
#'      )
#'      dem_params = validate_dem_params(dem_params, model_dims)
#'
#' }
#'
validate_dem_params <- function(dem_params, model_dims){
    required_names <- c("mort", "waa", "sexrat", "sel", "ret", "dmr", "surv_sel")
    required_dims <- c(model_dims$nyears, model_dims$nages, model_dims$nsexes, model_dims$nregions, model_dims$nfleets)

    input_names <- names(dem_params)
    missing_pars <- required_names[!(required_names %in% input_names)]
    if(length(missing_pars) > 0){
        stop(
            paste("Missing", length(missing_pars), "required demographic parameters:", paste0(missing_pars, collapse = ", "),".\nRequired params include:", paste0(required_names, collapse=", "), ".")
        )
    }

    for(n in input_names){
        d <- dem_params[[n]]
        dims <- dim(d)
        dim_match <- sum(dims != required_dims[1:length(dims)])
        if(dim_match > 0){
            stop(
                paste("Input matrix for", n, "is of wrong dimension.\nSupplied dimensions were", paste0(dims, collapse=", "), ", while required dimensions were", paste0(required_dims[1:length(dims)], collapse=", "), ".")
            )
        }
    }

    return(dem_params)
}
