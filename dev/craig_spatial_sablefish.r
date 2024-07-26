#' Recreating dynamics from SpatialSablefish Model_04
#' written and tested by Craig Marsh in 2022. Spatial
#' model that generated derived results can be found
#' at: https://craig44.github.io/SpatialSablefishAssessment/index.html.
#' 
#' All data and derived quantities used from the
#' SpatialSablefishAssessment model are preliminary.

rm(list=ls())
library(devtools)
devtools::load_all()

nyears <- length(1977:2020)
nages  <- 30
nsexes <- 2
nregions <- 5
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1977:2020,
    "age"  = 2:31,
    "sex"  = c("F", "M"),
    "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
    "fleet" = c("Fixed", "Trawl")
)

model_params <- set_model_params(nyears, nages, nsexes, nregions, nfleets)

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
load(file.path("data", "spatial_sablefish_assessment_data.rda"))

M <- spatial_sablefish_assessment_data$input_data$M[1,1] 
mort_matrix <- generate_param_matrix(M, dimension_names = dimension_names)

sexrat <- 0.5
sexrat_matrix <- generate_param_matrix(sexrat, dimension_names = dimension_names)

maturity <- array(t(spatial_sablefish_assessment_data$input_data$maturity), dim=c(nyears, nages), dimnames=dimension_names[c("time", "age")]) 
mat_matrix <- generate_param_matrix(maturity, dimension_names=dimension_names, by=c("time", "age"))

waa <- array(NA, dim=c(nyears, nages, nsexes), dimnames = dimension_names[c("time", "age", "sex")])
waa[,,1] <- t(spatial_sablefish_assessment_data$input_data$female_mean_weight_by_age)
waa[,,2] <- t(spatial_sablefish_assessment_data$input_data$male_mean_weight_by_age)
waa_matrix <- generate_param_matrix(waa, dimension_names = dimension_names, by=c("time", "age", "sex"))

sel <- array(NA, dim=c(nyears, nages, nsexes, nfleets), dimnames = dimension_names[c("time" ,"age", "sex", "fleet")])
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 0), , 1, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_f[,1]
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 0), , 2, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_m[,1]
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 1), , 1, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_f[,2]
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 1), , 2, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_m[,2]
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 2), , 1, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_f[,3]
sel[which(spatial_sablefish_assessment_data$input_data$fixed_sel_by_year_indicator == 2), , 2, 1] <- spatial_sablefish_assessment_data$outputs$sel_fixed_m[,3]
sel[,,1,2] <- t(spatial_sablefish_assessment_data$outputs$sel_trwl_f) 
sel[,,2,2] <- t(spatial_sablefish_assessment_data$outputs$sel_trwl_m)
sel_matrix <- generate_param_matrix(sel, dimension_names = dimension_names, by=c("time", "age", "sex", "fleet"))

ret <- 1
ret_matrix <- generate_param_matrix(ret, dimension_names = dimension_names, include_fleet_dim = TRUE)

dmr <- 0
dmr_matrix <- generate_param_matrix(dmr, dimension_names = dimension_names, include_fleet_dim = TRUE)

movement_matrix <- array(NA, dim=c(nregions, nregions, nages, nsexes), dimnames=dimension_names[c("region", "region", "age", "sex")])
movement_matrix[,,1:nages, 1:nsexes] <- spatial_sablefish_assessment_data$outputs$movement_matrix

dem_params <- list(
    waa=waa_matrix,
    mat=mat_matrix,
    mort=mort_matrix,
    sexrat=sexrat_matrix,
    sel=sel_matrix,
    ret=ret_matrix,
    dmr=dmr_matrix,
    movement=movement_matrix
)

#' 3. Generate initial numbers-at-age structure, recruitment timeseries,
#' and catch timeseries for the model.
init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames=c("time"=1977, dimension_names[c("age", "sex", "region")]))
init_naa[,,1,] <- spatial_sablefish_assessment_data$outputs$natage_f[,,1]
init_naa[,,2,] <- spatial_sablefish_assessment_data$outputs$natage_m[,,1]
recruitment_timeseries <- matrix(c(apply(spatial_sablefish_assessment_data$outputs$recruitment_yr, 1, sum), 0), nrow=nyears+1)

catch_timeseries <- apply(spatial_sablefish_assessment_data$input_data$fixed_fishery_catch+spatial_sablefish_assessment_data$input_data$trwl_fishery_catch, 2, sum)
catch_timeseries <- matrix(catch_timeseries, nrow=nyears)

#' 4. Apportion catch and recruitment across regions
rec_apportionment <- array(NA, dim=c(nyears+1, nregions), dimnames=list("time"=1977:2021, dimension_names[[c("region")]]))
rec_apportionment[1:(nyears),] <- t(apply(spatial_sablefish_assessment_data$outputs$recruitment_yr, 1, \(x) x/sum(x)))
rec_apportionment[nyears+1,] <- apply(rec_apportionment, 2, mean, na.rm=TRUE)

regional_catch <- spatial_sablefish_assessment_data$input_data$fixed_fishery_catch+spatial_sablefish_assessment_data$input_data$trwl_fishery_catch
catch_regional_apportionment <- t(apply(regional_catch, 2, \(x) x/sum(x)))

catch_fleet_apportionment <- array(NA, dim=c(nyears, nfleets, nregions), dimnames=dimension_names[c("time", "fleet", "region")])
catch_fleet_apportionment[,1,] <- t(spatial_sablefish_assessment_data$input_data$fixed_fishery_catch/regional_catch)
catch_fleet_apportionment[,2,] <- t(spatial_sablefish_assessment_data$input_data$trwl_fishery_catch/regional_catch)

#' 5. Setup Model Options
model_options <- setup_model_options(model_params)
model_options$region_apportionment <- catch_regional_apportionment
model_options$fleet_apportionment <- catch_fleet_apportionment
model_options$recruit_apportionment <- rec_apportionment
model_options$recruit_apportionment_random <- FALSE
model_options$do_recruits_move <- FALSE
model_options$simulate_observations <- FALSE

#' 6. Run operating model
om <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = recruitment_timeseries, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)


#' Plot spawning biomass
om_ssb <- compute_ssb(om$naa, dem_params)
dimnames(spatial_sablefish_assessment_data$outputs$SSB_yr) <- dimension_names[c('time', 'region')]

ssb_df <- reshape2::melt(om_ssb) %>% mutate(model="afscOM") %>%
    bind_rows(reshape2::melt(spatial_sablefish_assessment_data$outputs$SSB_yr) %>% mutate(model="craig")) %>%
    mutate(region = factor(
        region,
        levels = c("EGOA", "CGOA", "WGOA", "AI", "BS"),
        labels = c("East GOA", "Central GOA", "West GOA", "Aleutian Islands", "Bering Sea")
    ))

ggplot(ssb_df, aes(x=time, y=value, color=region, linetype=model))+
    geom_line(size=1)+
    scale_y_continuous(limits=c(0, 120))+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Spawning Biomass (mt)", title="Spawning Stock Biomass", color="Region", linetype="Model")+
    theme_bw()+
    facet_wrap(~region)


#' Plot recruitment
om_recruits <- apply(om$recruits[1:nyears,,,], c(1, 2), sum)
dimnames(om_recruits) <- dimension_names[c("time", "region")]
dimnames(spatial_sablefish_assessment_data$outputs$recruitment_yr) <- dimension_names[c("time", "region")]

rec_df <- reshape2::melt(om_recruits) %>% mutate(model="afscOM") %>%
    bind_rows(reshape2::melt(spatial_sablefish_assessment_data$outputs$recruitment_yr) %>% mutate(model="craig")) %>%
    mutate(region = factor(
        region,
        levels = c("EGOA", "CGOA", "WGOA", "AI", "BS"),
        labels = c("East GOA", "Central GOA", "West GOA", "Aleutian Islands", "Bering Sea")
    ))

ggplot(rec_df , aes(x=time, y=value, color=region, linetype=model))+
    geom_line(size=1)+
    scale_y_continuous(limits=c(0, 120))+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Recruits", title="Recruitment", color="Region", linetype="Model")+
    theme_bw()+
    facet_wrap(~region)


#' Plot Fs
om_f <- compute_total_f(om$faa)
dimnames(om_f) <- dimension_names[c("time", "region")]

spatial_assessment_faa <- array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
spatial_assessment_faa[,,1,,1] <- aperm(spatial_sablefish_assessment_data$outputs$F_fixed_f, perm=c(3, 1, 2))
spatial_assessment_faa[,,2,,1] <- aperm(spatial_sablefish_assessment_data$outputs$F_fixed_m, perm=c(3, 1, 2))
spatial_assessment_faa[,,1,,2] <- aperm(spatial_sablefish_assessment_data$outputs$F_trwl_f, perm=c(3, 1, 2))
spatial_assessment_faa[,,2,,2] <- aperm(spatial_sablefish_assessment_data$outputs$F_trwl_m, perm=c(3, 1, 2))
spatial_assessment_f <- compute_total_f(spatial_assessment_faa)

F_df <- reshape2::melt(om_f) %>% mutate(model="afscOM") %>%
    bind_rows(reshape2::melt(spatial_assessment_f) %>% mutate(model="craig")) %>%
    mutate(region = factor(
        region,
        levels = c("EGOA", "CGOA", "WGOA", "AI", "BS"),
        labels = c("East GOA", "Central GOA", "West GOA", "Aleutian Islands", "Bering Sea")
    ))

ggplot(F_df, aes(x=time, y=value, color=region, linetype=model))+
    geom_line(size=1)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="F", title="Fishing Mortality Rate", color="Region", linetype="Model")+
    theme_bw()+
    facet_wrap(~region)


#' Plot Catches
om_catches <- apply(om$land_caa, c(1, 4, 5), sum)
dimnames(om_catches) <- dimension_names[c("time", "region", "fleet")]

spatial_assessment_catches <- array(NA, dim=c(nyears, nregions, nfleets), dimnames=dimension_names[c("time", "region", "fleet")])
spatial_assessment_catches[,,1] <- t(spatial_sablefish_assessment_data$outputs$fixed_fishery_catch)
spatial_assessment_catches[,,2] <- t(spatial_sablefish_assessment_data$outputs$trwl_fishery_catch)

catch_df <- reshape2::melt(om_catches) %>% mutate(model="afscOM") %>%
    bind_rows(reshape2::melt(spatial_assessment_catches) %>% mutate(model="craig")) %>%
    mutate(region = factor(
        region,
        levels = c("EGOA", "CGOA", "WGOA", "AI", "BS"),
        labels = c("East GOA", "Central GOA", "West GOA", "Aleutian Islands", "Bering Sea")
    ))

ggplot(catch_df, aes(x=time, y=value, color=fleet, linetype=model))+
    geom_line(size=1)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Catch (mt)", title="Total Catch", color="Region", linetype="Model")+
    theme_bw()+
    facet_wrap(~region)


# Plot Age Comps
caas <- apply(om$caa[1:(nyears),,,,,drop=FALSE], c(1, 2, 4), sum)
dimnames(caas) <- dimension_names[c("time", "age", "region")]

reshape2::melt(caas) %>%
    as_tibble() %>%
    mutate(
        size_group = case_when(
                    age < 5 ~ "Small",
                    age < 9 ~ "Medium",
                    TRUE ~ "Large"
                )
    ) %>%
    group_by(region, size_group) %>%
    summarise(
        v = mean(sum(value))
    ) %>%
    tidyr::pivot_wider(names_from=size_group, values_from=v)
