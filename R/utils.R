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
#' @export
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
#' @export
#'
F_to_mu <- function(fy){
    return(1-exp(-fy))
}

#' Set variable names as names of list
#'
#' Automatically use variable names as list object names when
#' a sequence of variables is provided to a list constructor_
#'
#' @param ... R variables containing values to be put in a list
#'
#' @export
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
#' @param d dimension index, default: 1
#' @param drop controls whether singleton dimensions are dropped or preserved after subsetting default: TRUE
#' @return a list with the same elements as dem_params but containing
#' a single row from each list element
#'
#' @export
#'
#'
subset_dem_params <- function(dem_params, r, d=1, drop=TRUE){
    tmp <- rlang::duplicate(dem_params)
    ps <- names(tmp)
    for(n in ps){
        # if(any(is.na(dem_params[[n]]))) next;
        if(n == "movement"){
            if(d > 3) next; # Dont allow subsetting by region or shit goes to hell
            tmp[[n]] <- subset_matrix(dem_params[[n]], r, d+2, drop)
        }else{
            tmp[[n]] <- subset_matrix(dem_params[[n]], r, d, drop)
        }
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
#' @export
#'
#'
subset_matrix <- function(mat, r, d=1, drop=TRUE){
    if(is.null(mat)) return(NULL)
    if(is.vector(mat)) return(mat)
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
#' @export
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

        t <- do.call(afill, afill_params)
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
#' @export
#'
#' @examples
#' \dontrun{
#'      dimensions = list(nyears=60, nages=50, nsexes=2, nregions=1, nfleets=1)
#'      model_options = setup_model_options(dimensions)
#' }
#'
setup_model_options <- function(model_dimensions){

    return(
        list(
            removals_input = "catch",
            simulate_observations = TRUE,
            region_apportionment = matrix(1/model_dimensions$nregions, nrow=model_dimensions$nyears, ncol=model_dimensions$nregions),
            fleet_apportionment = array(1/model_dimensions$nfleets, dim=c(model_dimensions$nyears, model_dimensions$nfleets, model_dimensions$nregions)),
            recruit_apportionment = matrix(1/model_dimensions$nregions, nrow=(model_dimensions$nyears+1), ncol=model_dimensions$nregions),
            recruit_apportionment_random = FALSE,
            do_recruits_move = TRUE
        )
    )
}

#' Calculate Joint Selectivity and Retention Across Multiple Fleets
#'
#' Computes the average selectivity-at-age and retention-at-age acting
#' on a population when multiple fleets are present. Selectivity and
#' retention are weighted based on user supplied weights.
#'
#' @param sel selectivity-at-age ([1, nages, nesexes, nregions, nfleets])
#' @param ret retention-at-age ([1, nages, nsexes, nregions, nfleets])
#' @param prop_fs fleet weights
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_joint_selret(sel, ret, prop_fs=c(0.50, 0.50))
#' }
#'
calculate_joint_selret <- function(sel, ret, prop_fs=c(0.50, 0.50)){
    joint_self <- apply(sweep(sel[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_selm <- apply(sweep(sel[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_retf <- apply(sweep(ret[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_retm <- apply(sweep(ret[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))

    joint_sel <- array(NA, dim=dim(sel)[1:4])
    joint_sel[,,1,] <- joint_self
    joint_sel[,,2,] <- joint_selm

    joint_ret <- array(NA, dim=dim(ret)[1:4])
    joint_ret[,,1,] <- joint_retf
    joint_ret[,,2,] <- joint_retm

    return(list(sel=joint_sel, ret=joint_ret))
}


#' Generate output matrices
#' @param nyears number of years
#' @param nages number of ages
#' @param nsexes number of sexes
#' @param nregions number of regions
#' @param nfleets number of fishing fleets
#' @param nsurveys number of surveys
#' @export
generate_output_matrices <- function(nyears, nages, nsexes, nregions, nfleets, nsurveys){
    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
    # naa[1,,,] = init_naa

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))
    recruits    = array(NA, dim=c(nyears+1, 1, 1, nregions))

    survey_preds <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    survey_obs <- list(
        catch = array(NA, dim=c(nyears, 1, 1, nregions, nfleets)),
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    outputs <- listN(land_caa, disc_caa, caa, faa, naa, f, recruits, survey_preds, survey_obs)
    return(outputs)

}

#' update matrices
#'
#' @param output_matrices previous OM output
#' @param y year
#' @param input inputs (land_caa_tmp, caa_tmp, naa_tmp, etc.)
#' @param update_obs switch to also update the observation data, default: FALSE
#'
#' @export
#'
#' @examples
#' \dontrun{
#' update_output_matrices(output_matrices, y, input, update_obs=FALSE)
#' }
update_output_matrices <- function(output_matrices, y, input, update_obs=FALSE){

    # update state
    output_matrices$land_caa[y,,,,] <- input$land_caa_tmp
    output_matrices$disc_caa[y,,,,] <- input$disc_caa_tmp
    output_matrices$caa[y,,,,] <- input$caa_tmp
    output_matrices$faa[y,,,,] <- input$faa_tmp
    output_matrices$naa[y+1,,,] <- input$naa_tmp
    # output_matrices$recruits[y+1,,1,] <- apply(input$naa_tmp[,1,,,drop=FALSE], 3, sum)

    output_matrices$f[y,,,,] <- input$F_f_tmp
    # output_matrices$recruits[y+1,,,] <- rec

    if(update_obs){
        output_matrices$survey_preds$rpns[y,,,,] <- input$survey_preds$rpns
        output_matrices$survey_preds$rpws[y,,,,] <- input$survey_preds$rpws
        output_matrices$survey_preds$acs[y,,,,]  <- input$survey_preds$acs

        output_matrices$survey_obs$catch[y,,,,] <- input$survey_obs$catch
        output_matrices$survey_obs$rpns[y,,,,] <- input$survey_obs$rpns
        output_matrices$survey_obs$rpws[y,,,,] <- input$survey_obs$rpws
        output_matrices$survey_obs$acs[y,,,,]  <- input$survey_obs$acs
    }

    return(output_matrices)
}
