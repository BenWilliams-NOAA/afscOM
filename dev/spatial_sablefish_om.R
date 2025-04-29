rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(tictoc)
library(abind)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"

devtools::load_all(afscOM_dir)
# devtools::load_all(sablefishMSE_dir)

lapply(list.files("R", full.names = TRUE), source)

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons

spatial_sablefish_data <- readRDS("~/Desktop/Projects/SabieTMB/output/Spatial Assessment/data.RDS")
spatial_sablefish_estimates <- readRDS("~/Desktop/Projects/SabieTMB/output/Spatial Assessment/mle_rep.RDS")

spatial_sablefish_pars <- readRDS("~/Desktop/Projects/SabieTMB/output/Spatial Assessment/parameters.RDS")

load(file=file.path(here:::here(), "data", "spatial_sablefish_inputs.rda"))

nyears <- length(spatial_sablefish_inputs$years)
nages <- 30
nsexes <- 2
nregions <- 5
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:(nages+1),
    "sex"  = c("F", "M"),
    "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
    "fleet" = c("Fixed", "Trawl")
)

model_params <- afscOM::set_model_params(nyears, nages, nsexes, nregions, nfleets)
model_options <- afscOM::setup_model_options(model_params)

# 1. Generate demographic parameter matrices
mort_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_inputs$mort, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
waa_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_inputs$waa, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
mat_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_inputs$mat, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
sexrat_matrix <- afscOM::generate_param_matrix(0.5, dimension_names = dimension_names)

fish_selex_matrix = afscOM::generate_param_matrix(
    aperm(spatial_sablefish_inputs$fish_selex, perm=c(2, 3, 4, 1, 5)),
    by=c("time", "age", "sex", "region", "fleet"),
    dimension_names = dimension_names
)
fish_ret_matrix = afscOM::generate_param_matrix(1, dimension_names = dimension_names, include_fleet_dim = TRUE)
fish_dmr_matrix = afscOM::generate_param_matrix(0, dimension_names = dimension_names, include_fleet_dim = TRUE)
surv_selex_matrix = afscOM::generate_param_matrix(
    aperm(spatial_sablefish_inputs$surv_selex, perm=c(2, 3, 4, 1, 5)),
    by=c("time", "age", "sex", "region", "fleet"),
    dimension_names = c(dimension_names[1:4], list("fleet"=c("LL", "TRWL")))
)

movement <- spatial_sablefish_inputs$movement
dimnames(movement) <- c(rep(dimension_names[4], 2), dimension_names[1:3])

dem_params <- list(
    waa=waa_matrix,
    mat=mat_matrix,
    mort=mort_matrix,
    sexrat=sexrat_matrix,
    sel=fish_selex_matrix,
    ret=fish_ret_matrix,
    dmr=fish_dmr_matrix,
    surv_sel=surv_selex_matrix,
    movement=movement
)

# dp_y <- subset_dem_params(dem_params, 1, d=1, drop=FALSE)

# 2. Define initial population state
init_naa <- aperm(spatial_sablefish_inputs$naa[,1,,,drop=FALSE], perm=c(2, 3, 4, 1))
dimnames(init_naa) = c(list("time"=1), dimension_names[2:4])

# 3. Define recruitment timeseries
recruitment <- aperm(spatial_sablefish_inputs$recruitment, perm=c(2, 1))
dimnames(recruitment) <- dimension_names[c(1, 4)]

# 4. Define removals timeseries
fmort <- aperm(spatial_sablefish_inputs$fmort, perm=c(2, 3, 1))
catch <- aperm(spatial_sablefish_inputs$catch, perm=c(2, 3, 1))

# 5. Set Model Options
model_options$removals_input = "F"
model_options$simulate_observations = FALSE
model_options$do_recruits_move = FALSE

model_options$simulate_observations = TRUE
model_options$obs_pars <- list(
    # longline fishery, trawl fishery, longline survey, trawl survey
    is_survey   = c(0, 0, 1, 1),  # is this a survey (1) or fishery (0)
    qs          = c(1, 1, 6113.707, 9190.43), # catchability coefficient (q) for surveys
    catch_cv    = c(0.1, 0.1, 0, 0), # CV on catch for catch observations
    rpn         = c(0, 0, 1, 1), # should RPNs be computed (yes=1, no=0)
    rpn_cv      = c(0, 0, 0.1, 0.1), # RPN CV
    rpw         = c(0, 0, 1, 1), # should RPWs be computed (yes=1, no=0)
    rpw_cv      = c(0, 0, 0.1, 0.1), # RPW CV
    acs         = c(1, 1, 1, 1), # should age compositions be computed (yes=1, no=0)
    ac_samps    = c(200, 100, 200, 100), # total sample size for age composition observations
    ac_as_integers  = c(TRUE, TRUE, TRUE, TRUE), # return age comps as integers (TRUE) or proportions (FALSE)
    acs_agg_sex     = c(FALSE, FALSE, FALSE, FALSE) # should age comps be aggregated by sex
)

om_sim <- project(
    init_naa = init_naa, 
    removals_timeseries = fmort, 
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

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

matt_ssb <- aperm(spatial_sablefish_inputs$ssb, c(2, 1))
matt_catch <- apply(aperm(spatial_sablefish_inputs$catch, c(2, 3, 1)), c(1, 3), sum)
matt_bio <- apply(
    aperm(spatial_sablefish_inputs$naa[,1:nyears,,,drop=FALSE], c(2, 3, 4, 1))*dem_params$waa,
    c(1, 4),
    sum
)
matt_f <- apply(apply(aperm(spatial_sablefish_inputs$faa, c(2, 3, 4, 1, 5)), c(1, 4, 5), max), c(1, 2), sum)
dimnames(f) = list("time"=1:nyears, "region"=dimension_names[["region"]])
dimnames(matt_f) = dimnames(f)

plot_ssb(ssb, matt_ssb)
plot_catch(catch, matt_catch)
plot_bio(bio, matt_bio)
plot_f(f, matt_f)

plot_ssb(ssb)
plot_ssb(matt_ssb)
# spatial_sablefish_inputs <- list(
#     ssb = spatial_sablefish_estimates$SSB,
#     naa = spatial_sablefish_estimates$NAA,
#     caa = spatial_sablefish_estimates$CAA,
#     faa = spatial_sablefish_estimates$FAA,
#     fmort = spatial_sablefish_estimates$Fmort,
#     catch = spatial_sablefish_estimates$PredCatch,
#     recruitment = spatial_sablefish_estimates$Rec,
#     waa = spatial_sablefish_data$WAA,
#     mat = spatial_sablefish_data$MatAA,
#     mort = spatial_sablefish_estimates$natmort,
#     sexrat = spatial_sablefish_data$sexratio,
#     fish_selex = spatial_sablefish_estimates$fish_sel,
#     surv_selex = spatial_sablefish_estimates$srv_sel,
#     movement = spatial_sablefish_estimates$Movement,
#     years = spatial_sablefish_data$years
# )

# save(spatial_sablefish_inputs, file=file.path(here::here(), "data", "spatial_sablefish_inputs.rda"))

ll_rpn <- subset_matrix(om_sim$survey_preds$rpns, 1, d=5, drop=TRUE)
ak_ll_rpn <- apply(ll_rpn, 1, sum)

tw_rpn <- subset_matrix(om_sim$survey_preds$rpns, 2, d=5, drop=TRUE)
ak_tw_rpn <- apply(tw_rpn, 1, sum)

ll_ac_obs <- subset_matrix(om_sim$survey_obs$acs, 3, d=5, drop=TRUE)
ak_ll_ac_obs <- apply(ll_ac_obs, c(1, 2, 3), sum)

tw_rpw <- subset_matrix(om_sim$survey_preds$rpws, 2, d=5, drop=TRUE)
ak_tw_rpw <- apply(tw_rpw, c(1, 4), sum)

ak_rpw_props <- t(apply(ak_tw_rpw, 1, function(x) x/sum(x)))

dim(ll_ac_obs)
dim(ak_rpw_props)

test <- sweep(ll_ac_obs, c(1, 4), ak_rpw_props, FUN="*")
apply(round(apply(test, c(1, 2, 3), sum), 0), 1, sum)

ll_ac_obs
