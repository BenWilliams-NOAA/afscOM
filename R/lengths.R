compute_length_structure <- function(naa, sizeage_matrix){
    
    model_dims <- dim(sizeage_matrix)
    nlengths <- model_dims[1]
    nages <- model_dims[2]
    nsexes <- model_dims[3]
    nregions <- model_dims[4]
    nyears <- model_dims[5]
    
    v <- vapply(
        1:nyears,
        function(y){
            sizeage <- sizeage_matrix[,,,,y]
            sapply(
                1:model_params$nsexes,
                # Apply movement to the sexes individually
                function(s) naa[y,,s,] %*% sizeage_matrix[,,s,,y]
            )
        },
        FUN.VALUE = array(0, dim=c(nlengths, nsexes))
    )

    lengths <- array(aperm(v, c(3, 1, 2)), dim=c(nyears, nlengths, nsexes))
    return(lengths)
}
