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
#' @param include.fleet.dim whether to expand the output matrix to include a 5th dimension, indicating fleet structure
#'
#' @return a 4d or 5d array of dimensions (nyears, nages, nsexes, nregions, nfleets) filled across the specified dimensions by the specfied values.
#' 
#' @export generate_param_matrix
#'
#' @example
#'
generate_param_matrix <- function(vals, dimension.names, by=NA, include.fleet.dim=FALSE){

    # if(any(is.na(model.params)) & !exists("model.params")){
    #     stop("`model.params` object does not exist and was not provided.")
    # }

    # if(any(is.na(dimension.names)) & !exists("dimension.names")){
    #     stop("`dimension.names` object does not exist and was not provided.")
    # }

    model.params <- lapply(dimension.names, length)
    names(model.params) <- c("nyears", "nages", "nsexes", "nregions", "nfleets")[1:length(model.params)]

    if("fleet" %in% by){
        include.fleet.dim <- TRUE
    }

    dim.order <- c("time", "age", "sex", "region", "fleet")
    
    ndims <- 4
    if(include.fleet.dim){
        ndims <- 5
    }

    if(length(vals)==1){
        arr <- array(vals, dim=model.params[1:ndims], dimnames=dimension.names[1:ndims])
        return(arr)
    }

    tmp <- array(NA, dim=model.params[1:ndims], dimnames=dimension.names[1:ndims])
    val.dims <- dim(vals)

    # vals is a c vector
    if(is.null(val.dims)){
        val.dims <- length(vals)
    }
    n.dims <- val.dims
    expected.dims <- unlist(lapply(dimension.names[by], length))
    if(!all(n.dims == expected.dims)){
        stop(paste0("Wrong dimensions (", n.dims ,") for filling by ", by, ". Expected dimensions are ", expected.dims, ".\n"))
    }
    vals <- array(vals, dim=val.dims, dimnames=dimension.names[by])
    provided.dimensions <- dim.order %in% by
    afill.dimensions <- as.vector(!provided.dimensions, mode="list")
    for(i in which(provided.dimensions == TRUE)){
        afill.dimensions[[i]] <- missing_arg()
    }
    afill.params <- afill.dimensions[1:ndims]
    afill.params$x <- tmp
    afill.params$value <- vals

    tmp <- do.call("afill<-", afill.params)
    return(tmp)
}

dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
)

# Natural-Mortality-at-age
M <- 0.1
mort <- generate_param_matrix(M, dimension.names = dimension.names)

# Maturity-at-age
maa <- exp(-5.1560+0.7331*2:31)/(1+exp(-5.1560+0.7331*2:31))
maturity <- generate_param_matrix(maa, dimension.names = dimension.names, by="age")

# Weight-at-age
log.waa.f <- log(5.87)+3.02*log(1-exp(-0.17*(2:31+2.98)))
log.waa.m <- log(3.22)+3.02*log(1-exp(-0.27*(2:31+2.41)))
waa.f <- exp(log.waa.f)
waa.m <- exp(log.waa.m)

waa.matrix <- matrix(c(waa.f, waa.m), nrow=30)
colnames(waa.matrix) <- c("F", "M")
rownames(waa.matrix) <- 2:31

waa <- generate_param_matrix(waa.matrix, dimension.names = dimension.names, by=c("age", "sex"))

# Selectivity-at-age per fleet
all_selex <- dget("../sablefish-mse/scratch/tem_2023.rdat")$agesel
selex.f.ll_fish <- all_selex[,"fish1sel.f"]
selex.m.ll_fish <- all_selex[,"fish1sel.m"]
selex.f.tw_fish <- all_selex[,"fish3sel.f"]
selex.m.tw_fish <- all_selex[,"fish3sel.m"]

selex.arr <- array(NA, dim=c(30, 2, 2), dimnames=list(2:31, c("F", "M"), c("fixed", "trawl")))
selex.arr[,1,1] <- selex.f.ll_fish
selex.arr[,2,1] <- selex.m.ll_fish
selex.arr[,1,2] <- selex.f.tw_fish
selex.arr[,2,2] <- selex.m.tw_fish

sel <- generate_param_matrix(selex.arr, dimension.names = dimension.names, by=c("age", "sex", "fleet"))

l <- list(
    single = mort,
    age = maturity,
    age.sex = waa,
    age.sex.fleet = sel
)

saveRDS(l, "data/good_dem_matrices.RDS")
