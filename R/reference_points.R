compute_naapr <- function(nages, mort, mat, waa, sel, ret, F){
    naa <- rep(NA, nages)
    naa[1] <- 1
    zaa <- mort + sel*ret*F
    for(a in 2:(nages)){
        naa[a] <- naa[a-1]*exp(-zaa[a-1])
    }
    return(naa)
}

#' Compute Spawning Biomass per Recruit (SBPR)
#'
#' Compute SBPR under a given level of fishing mortality.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantenous fishing mortality rate
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compute_sbpr(nages, mort, mat, waa, sel, ret, F)
#' }
#'
#'
compute_sbpr <- function(nages, mort, mat, waa, sel, ret, F){
    naa <- compute_naapr(nages, mort, mat, waa, sel, ret, F)
    ssb <- sum(naa*mat*waa, na.rm = TRUE)
    return(ssb)
}

#' Compute Spawning Potential Ratio (SPR)
#'
#' Compute SPR, the ratio of spawning biomass under a
#' given leve of fishing mortality relative to unfished
#' spawning biomass.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantenous fishing mortality rate
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compute_spr(nages, mort, mat, waa, sel, ret, F)
#' }
compute_spr <- function(nages, mort, mat, waa, sel, ret, F){
    ssb_unfished <- compute_sbpr(nages, mort, mat, waa, sel, ret, F=0)
    ssb_fished   <- compute_sbpr(nages, mort, mat, waa, sel, ret, F)
    return(ssb_fished/ssb_unfished)
}

#' Find F that yields a given SPR%
#'
#' Use bisection algorithm to identify the level of
#' fishing mortality required to yield an SPR of x%.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param target_x desired SPR proportion
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compute_fx(nages, mort, mat, waa, sel, ret, target_x=0.35)
#' }
compute_fx <- function(nages, mort, mat, waa, sel, ret, target_x=0.35){
    range <- c(0,2)
    n.iter <- 20
    i <- 1
    for(i in 1:n.iter) {
      midpoint <- mean(range)
      spr <- compute_spr(nages, mort, mat, waa, sel, ret, F=midpoint)
      if(spr > target_x) {
        range[1] <- midpoint
        }else {
        range[2] <- midpoint
      }
    }
    Fx <- midpoint
    return(Fx)
}

#' Compute average SSB under a level of F
#'
#' Compute average SSB under a given level of
#' fishing mortality.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantaneous fishing mortality rate
#' @param avg_rec average recruitment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compute_bx(nages, mort, mat, waa, sel, ret, F, avg_rec)
#' }
compute_bx <- function(nages, mort, mat, waa, sel, ret, F, avg_rec){
  return(avg_rec*compute_sbpr(nages, mort, mat, waa, sel, ret, F))
}

#' Calculate NPFMC groundfish reference points
#'
#' Calculate F_OFL (F_35%), F_ABC (F_40%), B40
#' (SSB under F_ABC), and B0 (unfished SSB) reference
#' points, using the traditional AFSC approach.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param avg_rec average recruitment
#'
#' @return list of F40, F35, and B40
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_npfmc_ref_points(nages, mort, mat, waa, sel, ret, avg_rec)
#' }
calculate_npfmc_ref_points <- function(nages, mort, mat, waa, sel, ret, avg_rec){
  F35 <- compute_fx(nages, mort, mat, waa, sel, ret, target_x=0.35)
  F40 <- compute_fx(nages, mort, mat, waa, sel, ret, target_x=0.40)
  B40 <- compute_bx(nages, mort, mat, waa, sel, ret, F=F40, avg_rec=avg_rec)
  B0 <- compute_bx(nages, mort, mat, waa, sel, ret, F=0, avg_rec=avg_rec)
  return(list(F40=F40, F35=F35, B40=B40, B0))
}


#' Calculate F/B Reference Points via SPR Methods
#'
#' Calculate fishing mortality and biological
#' fisheries reference points via traditional SPR
#' methods.
#'
#' @param nages number of ages in age structure
#' @param dem_params list of demographic parameters
#' subset to a single year
#' @param joint_selret list of total selectivity and
#' retention across all fleets, weighted by relative
#' fleet-specific fishing mortality (compute using
#' `calculate_joint_selret`)
#' @param spr_target desired SPR proportion
#' @param rec recruitment used to calculate BRP
#'
#' @return list of Fref and Bref
#' @export
#'
#' @examples
#' \dontrun{
#' calculate_spr_refpoints(nages, dem_params, joint_selret, spr_target, rec)
#' }
calculate_spr_refpoints <- function(nages, dem_params, joint_selret, spr_target, rec){
    mort <- dem_params$mort[,,1,]
    mat <- dem_params$mat[,,1,]
    waa <- dem_params$waa[,,1,]
    sel <- joint_selret$sel[,,1,,drop=FALSE]
    ret <- joint_selret$ret[,,1,,drop=FALSE]

    Fref <- compute_fx(nages, mort, mat, waa, sel, ret, target_x = spr_target)
    Bref <- compute_bx(nages, mort, mat, waa, sel, ret, F=Fref, avg_rec=rec)
    return(listN(Fref, Bref))
}
