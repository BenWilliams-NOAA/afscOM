#' Setup model parameters
#'
#' @param nyears years of input data
#' @param nages number of ages
#' @param nsexes number of sexes
#' @param nregions number of regions
#' @param nfleets number of fishing fleets
#'
#' @export
#'
#' @examples
#' /dontrun{
#' set_model_params(nyears, nages, nsexes=1, nregions=1, nfleets=1)
#' }
set_model_params <- function(nyears, nages, nlengths=1, nsexes=1, nregions=1, nfleets=1){
    return(
        list(
            nyears = nyears,
            nages = nages,
            nlengths = nlengths,
            nsexes = nsexes,
            nregions = nregions,
            nfleets = nfleets
        )
    )
}
