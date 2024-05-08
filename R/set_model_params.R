#' Setup model parameters
#'
#' @param nyears years of input data
#' @param nages number of ages
#' @param nsexes number of sexes
#' @param nregions number of regions
#' @param nfleets number of fishing fleets
#'
#' @return
#' @export
#'
#' @examples
#' /dontrun{
#' set_model_params(nyears, nages, nsexes=1, nregions=1, nfleets=1)
#' }
set_model_params <- function(nyears, nages, nsexes=1, nregions=1, nfleets=1){
    return(
        list(
            nyears = nyears,
            nages = nages,
            nsexes = nsexes,
            nregions = nregions,
            nfleets = nfleets
        )
    )
}
