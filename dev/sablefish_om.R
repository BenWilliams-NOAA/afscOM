rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
library(afscOM)
#library(devtools)
#devtools::load_all()

assessment <- dget("data/test.rdat")

#' 1. Define model dimensions and dimension names
#'  - nyears: number of years over which to simulate
#'  - nages: number of age classes to use
#'  - nsexes: number of sexes (1 = females only)
#'  - nregions: number of spatial regions
#'  - nfleets: number of fishing fleets (survey fleets are handled seperately)
#'
#' Dimension names must be defined in order to use the helper function
#' `generate_param_matrix` which handles filling complex multi-dimensional
#' arrays with smaller dimensions values, vectors, or matrices.
nyears <- 64
nages  <- 30
nsexes <- 2
nregions <- 1
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:31,
    "sex"  = c("F", "M"),
    "region" = "alaska",
    "fleet" = c("Fixed", "Trawl")
)

model_params <- set_model_params(nyears, nages, nsexes, nregions, nfleets)

#' 2. Generate list with appropiate demographic parameters
#' Required parameters are:
#'  - waa: weight-at-age
#'  - mat: maturity-at-age
#'  - mort: natural-mortality-at-age
#'  - sexrat: sex ratio
#'  - sel: selectivity-at-age
#'  - ret: retention-at-age
#'  - dmr: discard-mortality-at-age
#'  - surv_sel: survey selectivity-at-age

M <- 0.113179
mort <- generate_param_matrix(M, dimension_names = dimension_names)

prop_males <- 0.5
sexrat <- generate_param_matrix(prop_males, dimension_names = dimension_names)

retention <- 1.0
ret <- generate_param_matrix(retention, dimension_names = dimension_names, include_fleet_dim = TRUE)

discard <- 0.0
dmr <- generate_param_matrix(discard, dimension_names = dimension_names, include_fleet_dim = TRUE)

maturity <- assessment$growthmat[, "mage.block1"]
mat <- generate_param_matrix(maturity, dimension_names = dimension_names, by="age")

weight_mat <- matrix(NA, nrow=nages, ncol=nsexes, dimnames=dimension_names[c("age", "sex")])
weight_mat[,1] <- assessment$growthmat[, "wt.f.block1"]
weight_mat[,2] <- assessment$growthmat[, "wt.m.block1"]
waa <- generate_param_matrix(weight_mat, dimension_names = dimension_names, by=c("age", "sex"))

selex_mat <- array(NA, dim=c(nages, nsexes, nfleets), dimnames=dimension_names[c("age", "sex", "fleet")])
selex_mat[,1,1] <- assessment$agesel[,"fish1sel.f"]
selex_mat[,2,1] <- assessment$agesel[,"fish1sel.m"]
selex_mat[,1,2] <- assessment$agesel[,"fish3sel.f"]
selex_mat[,2,2] <- assessment$agesel[,"fish3sel.m"]
sel <- generate_param_matrix(selex_mat, dimension_names = dimension_names, by=c("age", "sex", "fleet"), include_fleet_dim = TRUE)
sel[(36:56),,1,,1] <- matrix(rep(assessment$agesel[, "fish4sel.f"], length(36:56)), ncol=nages, byrow=TRUE)
sel[(36:56),,2,,1] <- matrix(rep(assessment$agesel[, "fish4sel.m"], length(36:56)), ncol=nages, byrow=TRUE)
sel[(57:64),,1,,1] <- matrix(rep(assessment$agesel[, "fish5sel.f"], length(57:64)), ncol=nages, byrow=TRUE)
sel[(57:64),,2,,1] <- matrix(rep(assessment$agesel[, "fish5sel.m"], length(57:64)), ncol=nages, byrow=TRUE)

survey_selex_mat <- array(NA, dim=c(nages, nsexes, 2), dimnames=dimension_names[c("age", "sex", "fleet")])
survey_selex_mat[,1,1] <- assessment$agesel[,"srv1sel.f"]
survey_selex_mat[,2,1] <- assessment$agesel[,"srv1sel.m"]
survey_selex_mat[,1,2] <- assessment$agesel[,"srv7sel.f"]
survey_selex_mat[,2,2] <- assessment$agesel[,"srv7sel.m"]
survey_sel <- generate_param_matrix(survey_selex_mat, dimension_names = dimension_names, by=c("age", "sex", "fleet"), include_fleet_dim = TRUE)
survey_sel[(57:64),,1,,1] <- matrix(rep(assessment$agesel[, "srv10sel.f"], length(57:64)), ncol=nages, byrow=TRUE)
survey_sel[(57:64),,2,,1] <- matrix(rep(assessment$agesel[, "srv10sel.m"], length(57:64)), ncol=nages, byrow=TRUE)

dem_params <- list(
    waa=waa,
    mat=mat,
    mort=mort,
    sexrat=sexrat,
    sel=sel,
    ret=ret,
    dmr=dmr,
    surv_sel=survey_sel
)

#dem_params <- validate_dem_params(dem_params, model_params)

#' 3. Define starting population condition (NAA)
#' Here, we will use the NAA in 1960 from the 2023
#' Alaska sablefish assessment.
init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames = list(1, 2:31, c("F", "M"), "alaska"))
init_naa[,,1,] <- assessment$natage.female["1960",]
init_naa[,,2,] <- assessment$natage.male["1960",]

#' 4. Define recruitment timeseries
#' Where possible (e.g. where recruitment is independent of
#' the population state), recruitment should be defined
#' external to the OM simulation loop to ensure reproducability
#' across model runs.
#'
#' Here we using the recruitment timeseries from the 2023 Alaska
#' sablefish assessment.

# Multiply by 2 because two sexes
recruitment <- assessment$natage.female[,1]*2
recruitment <- c(recruitment, recruitment[64])


#' 5. Define catch history
#' Catch histories are defined by three variables:
#'  - TAC: the total allowable catch each year across the entire model
#'  - region apportionment: the proportion of the TAC allocated to each
#'                          region in the model
#'  - fleet apportionment: the proportion of the TAC allocation to each
#'                         fishing fleet in each region in the model
#' The apportionment timeseries are placed in a `model_options` list.

TACs <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
fixed_fleet_prop <- assessment$t.series[,"Catch_HAL"]/(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
trawl_fleet_prop <- 1-fixed_fleet_prop

model_options <- list(
    region_apportionment = list(1),
    fleet_apportionment = list(fixed_fleet_prop, trawl_fleet_prop)
)

#' 6. Define parameters for observation processes
#' Observation process parameters include catchability coefficients 
#' (q), observation errors, and sample sizes for age/length comps.

obs_pars <- list(
    surv_ll = list(
        q = 6.41538,
        rpn_cv = 0.20,
        rpw_cv = 0.10,
        ac_samps = 1000
    ),
    surv_tw = list(
        q = 0.8580,
        rpw_cv = 0.10,
        ac_samps = 1000
    ),
    fish_fx = list(
        ac_samps = 1000
    )
)

model_options$obs_pars <- obs_pars

#' 6. Setup empty array to collect derived quantities from the OM
#' It is left to the user to decide what information to store, and
#' in what format they would like to store it.
#'
#' Here, we are storing all OM outputs as they are returned from the
#' OM.
land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
naa[1,,,] = init_naa

survey_preds <- list(
    ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_ac = array(NA, dim=c(nyears, nages, nsexes, nregions)),
    fxfish_caa = array(NA, dim=c(nyears, nages, nsexes, nregions))
)

survey_obs <- list(
    ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions)),
    fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions))
)

#' 7. Run the OM forward in time
#' The `project` function handles projecting the population
#' forward in time 1 timestep based on the provded TAC, current
#' population size, provided recruitment, and other model
#' options. In order to simulate multiple years, the `project`
#' function MUST be wrapped in a loop over the intended number
#' of simulation years.
#'
#' The `project` function does not "know" about the presence of
#' an external "year loop", and thus expects arguments to
#' represent data for a single year at a time. Helper functions
#' have been provided to help with correctly subsetting the
#' demographic matrices to ensure input data is of the correct
#' dimensionality.

for(y in 1:nyears){

    # Subset the demographic parameters list to only the current year
    # and DO NOT drop lost dimensions.
    dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
    fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))
    out_vars <- project(
        TAC=TACs[y],
        dem_params=dp.y,
        prev_naa=naa[y,,,, drop = FALSE],
        recruitment=recruitment[y+1],
        fleet.props = fleet.props,
        options=model_options
    )

    # update state
    land_caa[y,,,,] <- out_vars$land_caa_tmp
    disc_caa[y,,,,] <- out_vars$disc_caa_tmp
    caa[y,,,,] <- out_vars$caa_tmp
    faa[y,,,,] <- out_vars$faa_tmp
    naa[y+1,,,] <- out_vars$naa_tmp

    survey_preds$ll_rpn[y,,,] <- out_vars$surv_preds$ll_rpn
    survey_preds$ll_rpw[y,,,] <- out_vars$surv_preds$ll_rpw
    survey_preds$tw_rpw[y,,,] <- out_vars$surv_preds$tw_rpw
    survey_preds$ll_ac[y,,,] <- out_vars$surv_preds$ll_ac
    survey_preds$fxfish_caa[y,,,] <- out_vars$surv_preds$fxfish_caa

    survey_obs$ll_rpn[y,,,] <- out_vars$surv_obs$ll_rpn
    survey_obs$ll_rpw[y,,,] <- out_vars$surv_obs$ll_rpw
    survey_obs$tw_rpw[y,,,] <- out_vars$surv_obs$tw_rpw
    survey_obs$ll_acs[y,,,] <- out_vars$surv_obs$ll_ac_obs
    survey_obs$fxfish_acs[y,,,] <- out_vars$surv_obs$fxfish_caa_obs

}

#' 8. Plot OM Results
#'

ssb <- apply(naa[1:64,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum)
bio <- apply(naa[1:64,,,]*dem_params$waa[,,,], 1, sum)
catch <- apply(caa, 1, sum)
f <- apply(apply(faa, c(1, 5), \(x) max(x)), 1, sum)

ssb_comp <- data.frame(
    year=1960:2023,
    assess_ssb=assessment$t.series[, "spbiom"],
    om_ssb=ssb,
    assess_catch=TACs,
    om_catch=catch,
    assess_bio=assessment$t.series[, "totbiom"],
    om_bio = bio,
    assess_f = assessment$t.series[,"fmort"],
    om_f = f
)

library(ggplot2)

p1 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=assess_ssb), col="black", size=0.7)+
    geom_line(aes(y=om_ssb), col="red", size=0.7)+
    scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    labs(y="SSB", x="Year", title="Spawning Biomass Comparison")+
    theme_bw()

p2 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=assess_catch), col="black", size=0.7)+
    geom_line(aes(y=om_catch), col="red", size=0.7)+
    scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    labs(y="Catch", x="Year", title="Total Catch Comparison")+
    theme_bw()

p3 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=assess_bio), col="black", size=0.7)+
    geom_line(aes(y=om_bio), col="red", size=0.7)+
    scale_y_continuous(limits=c(0, 750), breaks=seq(0, 750, 100))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    labs(y="Biomass", x="Year", title="Total Biomass Comparison")+
    theme_bw()

p4 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=assess_f), col="black", size=0.7)+
    geom_line(aes(y=om_f), col="red", size=0.7)+
    scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.2, 0.05))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    labs(y="F", x="Year", title="Fishing Mortality Comparison")+
    theme_bw()

library(gridExtra)

p <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)
ggsave("~/Desktop/sablefish_assess_om.png", plot=p, width=8, height=8, units=c("in"))

(faa[,,1,1,1] - assessment$faa.fish1.f)/assessment$faa.fish1.f
# caa[,,1,1,1] - assessment$
y=30
faa[y,,1,1,1] - assessment$faa.fish1.f[y,]

100*(ssb - assessment$t.series[,"spbiom"])/assessment$t.series[,"spbiom"]

 - assessment$t.series[,"fmort"]
