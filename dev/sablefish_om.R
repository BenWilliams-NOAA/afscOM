rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
#library(afscOM)
library(devtools)
devtools::load_all()

load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
assessment <- sablefish_assessment_data

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
model_options <- setup_model_options(model_params)

#' 2. Generate list with appropriate demographic parameters
#' Required parameters are:
#'  - waa: weight-at-age
#'  - mat: maturity-at-age
#'  - mort: natural-mortality-at-age (instantaneous)
#'  - sexrat: sex ratio
#'  - sel: selectivity-at-age
#'  - ret: retention-at-age
#'  - dmr: discard-mortality-at-age (instantaneous)
#'  - surv_sel: survey selectivity-at-age

M <- 0.113179
mort <- generate_param_matrix(M, dimension_names = dimension_names)
mort[,,2,] <- mort[,,2,]-0.00819813

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

# saveRDS(dem_params, file="data/sabelfish_dem_params.RDS")

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
#'
#' Alternatively, a timeseries of fishing mortality rates, F, can be
#' provided in lieu of TACs. In this case, the timeseries of F for each
#' fleet, in each area, must be provided in an array of dimensions
#' [nyears, 1, 1, nregions, nfleets].
#'
#' Always specify the `removals_input` in the `model_options` list.
#' `removals_input` can be either "catch" or "F", depending on what
#' data you are providing as removals input.

TACs <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
# fixed_fleet_prop <- assessment$t.series[,"Catch_HAL"]/(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
# trawl_fleet_prop <- 1-fixed_fleet_prop

model_options <- setup_model_options(model_dimensions = model_params)
model_options$random_apportion_recruits <- FALSE
# model_options <- list(
#     region_apportionment = matrix(1, nrow=nyears, ncol=nregions),
#     fleet_apportionment = matrix(rep(c(fixed_fleet_prop, trawl_fleet_prop), nyears), ncol=nfleets),
#     removals_input = "catch"
# )


f_timeseries <- assessment$t.series[,c("F_HAL", "F_TWL")] %>% as.matrix
f_timeseries <- array(f_timeseries, dim=c(nyears, 1, 1, 1, 2),
                dimnames = list("time"=1:nyears,
                                age="all",
                                sex="all",
                                "region"="alaska",
                                "fleet"=c("Fixed", "Trawl")))

model_options$removals_input = "F"

#' 6. Define parameters for observation processes
#' Observation process parameters include catchability coefficients
#' (q), observation errors, and sample sizes for age/length comps.

obs_pars <- list(
    # longline fishery, trawl fishery, longline survey, trawl survey
    is_survey   = c(0, 0, 1, 1),  # is this a survey (1) or fishery (0)
    qs          = c(1, 1, 6.41, 0.85), # catchability coefficient (q) for surveys
    catch_cv    = c(0.1, 0.1, 0, 0), # CV on catch for catch observations
    rpn         = c(0, 0, 1, 1), # should RPNs be computed (yes=1, no=0)
    rpn_cv      = c(0, 0, 0.1, 0.1), # RPN CV
    rpw         = c(0, 0, 1, 1), # should RPWs be computed (yes=1, no=0)
    rpw_cv      = c(0, 0, 0.1, 0.1), # RPW CV
    acs         = c(1, 1, 1, 1), # should age compositions be computed (yes=1, no=0)
    ac_samps    = c(50, 30, 50, 30), # total sample size for age composition observations
    ac_as_integers  = c(TRUE, TRUE, TRUE, TRUE), # return age comps as integers (TRUE) or proportions (FALSE)
    acs_agg_sex     = c(FALSE, FALSE, FALSE, FALSE) # should age comps be aggregated by sex
)
model_options$simulate_observations <- TRUE
model_options$obs_pars <- obs_pars

#' 6. Run the OM forward in time
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
om_sim <- project_multi(
            init_naa = init_naa, 
            removals_timeseries = f_timeseries, 
            recruitment = recruitment, 
            dem_params = dem_params, 
            nyears = nyears, 
            model_options = model_options
          )

#' 8. Plot OM Results
#'

dimnames(om_sim$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="alaska")
dimnames(om_sim$caa) <- list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"=c("Fixed", "Trawl"))

# ssb <- apply(om_sim$naa[1:64,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum)
# bio <- apply(om_sim$naa[1:64,,,]*dem_params$waa[,,,], 1, sum)
# catch <- apply(om_sim$caa, 1, sum)
# f <- apply(apply(om_sim$faa, c(1, 5), \(x) max(x)), 1, sum)
ssb <- compute_ssb(om_sim$naa, dem_params)
bio <- compute_bio(om_sim$naa, dem_params)
catch <- compute_total_catch(om_sim$caa)
fleet_catch <- compute_fleet_catch(om_sim$caa)
f <- compute_total_f(om_sim$faa)

library(patchwork)

# Plot SSB, catch, biomass, and fishing mortality
p1 <- plot_ssb(ssb)
p2 <- plot_catch(catch)
p3 <- plot_bio(bio)
p4 <- plot_f(f)
(p1+p2)/(p3+p4)+plot_layout(guides="collect")

# Plot NAA and CAA
p5 <- plot_atage(om_sim$naa) + labs(title="Numbers-at-age")
p6 <- plot_atage(om_sim$caa) + labs(title="Catch-at-age")
(p5+p6) + plot_layout(guides="collect")

# Calculate and plot timeseries of B40
f_props <- aperm(apply(apply(om_sim$faa, c(1, 5), \(x) max(x)), 1, \(y) y/sum(y)), c(2, 1))
brps <- sapply(
    1:nyears,
    function(x){
        dp_y <- subset_dem_params(dem_params, r=x, d=1, drop=FALSE)
        selret <- calculate_joint_selret(dp_y$sel, dp_y$ret, prop_fs = f_props[x,])
        rec <- unlist(slider::slide(om_sim$naa[,1,1,], mean, .before=Inf))
        rps <- calculate_spr_refpoints(nages, dp_y, selret, spr_target=0.40, rec[x])
        return(rps$Bref)
    }
)


# Plot survey index observations
ll_surv_data <- data.frame(assessment$obssrv3) %>% rownames_to_column("Year")
ll_surv_data$Year <- as.numeric(ll_surv_data$Year)
ll_surv_data$om_pred <- om_sim$survey_preds$rpn[31:64,,,,1]
ll_surv_data$om <- om_sim$survey_obs$rpn[31:64,,,,1]
ll_surv_data$om.lci <- ll_surv_data$om -1.96*0.20*ll_surv_data$om
ll_surv_data$om.uci <- ll_surv_data$om +1.96*0.20*ll_surv_data$om

p1 <- ggplot(ll_surv_data, aes(x=Year, y=obssrv3, group=1))+
    geom_pointrange(aes(ymin=obssrv3.lci, ymax=obssrv3.uci, color="Assessment"))+
    geom_line(aes(y=predsrv3, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 3000), breaks=seq(0, 3000, 500))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="LL Survey RPN", x="Year", title="LL Survey RPN Comparison")+
    theme_bw()

ll_surv_data <- data.frame(assessment$obssrv1) %>% rownames_to_column("Year")
ll_surv_data$Year <- as.numeric(ll_surv_data$Year)
ll_surv_data$om_pred <- om_sim$survey_preds$rpw[31:64,,,,1]
ll_surv_data$om <- om_sim$survey_obs$rpw[31:64,,,,1]
ll_surv_data$om.lci <- ll_surv_data$om -1.96*0.10*ll_surv_data$om
ll_surv_data$om.uci <- ll_surv_data$om +1.96*0.10*ll_surv_data$om

p2 <- ggplot(ll_surv_data, aes(x=Year, y=obssrv1, group=1))+
    geom_pointrange(aes(ymin=obssrv1.lci, ymax=obssrv1.uci, color="Assessment"))+
    geom_line(aes(y=predsrv1, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 5500), breaks=seq(0, 5500, 1000))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="LL Survey RPW", x="Year", title="LL Survey RPW Comparison")+
    theme_bw()

tw_surv_data <- data.frame(assessment$obssrv7) %>% rownames_to_column("Year")
tw_surv_data$Year <- as.numeric(tw_surv_data$Year)
tw_surv_data$om_pred <- om_sim$survey_preds$rpw[tw_surv_data$Year-1960+1,,,,2]
tw_surv_data$om <- om_sim$survey_obs$rpw[tw_surv_data$Year-1960+1,,,,2]
tw_surv_data$om.lci <- tw_surv_data$om -1.96*0.10*tw_surv_data$om
tw_surv_data$om.uci <- tw_surv_data$om +1.96*0.10*tw_surv_data$om

p3 <- ggplot(tw_surv_data, aes(x=Year, y=obssrv7, group=1))+
    geom_pointrange(aes(ymin=obssrv7.lci, ymax=obssrv7.uci, color="Assessment"))+
    geom_line(aes(y=predsrv7, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 500), breaks=seq(0, 500, 100))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="TW Survey RPW", x="Year", title="TW Survey RPW Comparison")+
    theme_bw()

p1/p2/p3 + plot_layout(guides="collect")


survey_obs$ll_rpn
survey_preds$ll_rpn
n = length(survey_obs$ll_rpn)

exp((1/n)*sum(log(survey_obs$ll_rpn/pn[1:64])))
ssb_comp$om_bio

pn <- apply(naa[1:64,,,]*dem_params$surv_sel[1:64,,,,1], 1, sum)
