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
