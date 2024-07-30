rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
#library(afscOM)
library(devtools)
devtools::load_all()

load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
assessment <- sablefish_assessment_data

ssb_recuit_apportionment <- function(naa, dem_params){
    ssb <- naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,]
    ssb_by_region <- apply(ssb, 2, sum)
    return(
        array(ssb_by_region/sum(ssb_by_region), dim=c(1, ncol(ssb)))
    )
}


age2_apportionment <- function(naa, dem_params){
    age2 <- naa[,1,,]
    age2_by_region <- apply(age2, 2, sum)
    return(
        array(age2_by_region/sum(age2_by_region), dim=c(1, ncol(age2)))
    )
}

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
nyears <- 500
nages  <- 30
nsexes <- 2
nregions <- 5
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:31,
    "sex"  = c("F", "M"),
    "region" = c("EGOA", "CGOA", "WGOA", "BS", "AI"),
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

# Movement matrices from Hanselman et al. 2018
small_movement <- matrix(
    c(
        0.503, 0.294, 0.127, 0.021, 0.019, 
        0.372, 0.325, 0.180, 0.057, 0.053, 
        0.271, 0.304, 0.196, 0.112, 0.110,
        0.070, 0.148, 0.172, 0.567, 0.042,
        0.038, 0.085, 0.105, 0.049, 0.722
    ),
    nrow=nregions, ncol=nregions, byrow=TRUE
)
small_movement <- t(apply(small_movement, 1, \(x) x/sum(x)))

medium_movement <- matrix(
    c(
        0.584, 0.261, 0.076, 0.014, 0.021, 
        0.369, 0.356, 0.139, 0.049, 0.075, 
        0.271, 0.339, 0.151, 0.091, 0.140, 
        0.081, 0.200, 0.151, 0.502, 0.065, 
        0.073, 0.183, 0.141, 0.054, 0.548
    ),
    nrow=nregions, ncol=nregions, byrow=TRUE
)
medium_movement <- t(apply(medium_movement, 2, \(x) x/sum(x)))

large_movement <- matrix(
    c(
       0.550, 0.272, 0.094, 0.023, 0.024, 
       0.458, 0.306, 0.114, 0.050, 0.055, 
       0.423, 0.304, 0.117, 0.067, 0.076, 
       0.172, 0.227, 0.115, 0.395, 0.087, 
       0.153, 0.207, 0.106, 0.030, 0.501 
    ),
    nrow=nregions, ncol=nregions, byrow=TRUE
)
large_movement <- t(apply(large_movement, 2, \(x) x/sum(x)))

movement_matrix <- array(NA, dim=c(nregions, nregions, nages, nsexes))
movement_matrix[,,1:3,] <- small_movement
movement_matrix[,,4:6,] <- medium_movement
movement_matrix[,,7:nages,] <- large_movement

dem_params <- list(
    waa=waa,
    mat=mat,
    mort=mort,
    sexrat=sexrat,
    sel=sel,
    ret=ret,
    dmr=dmr,
    surv_sel=survey_sel,
    movement=movement_matrix
)


init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames = list(1, 2:31, c("F", "M"), c("EGOA", "CGOA", "WGOA", "BS", "AI")))
init_naa[,,1,] <- assessment$natage.female["1960",]/nregions
init_naa[,,2,] <- assessment$natage.male["1960",]/nregions

recruitment <- assessment$natage.female[,1]*2
recruitment <- as.vector(sample(recruitment, size=nyears+1, replace=TRUE))

catch_timeseries <- matrix(rep(0, nyears), nrow=nyears)

model_options <- list()
model_options$removals_input = "catch"
model_options$region_apportionment = matrix(rep(rep(1/6, nregions), nyears), ncol=nregions)
model_options$fleet_apportionment = matrix(1, nrow=nyears, ncol=nfleets)
obs_pars <- list(
    # longline fishery, trawl fishery, longline survey, trawl survey
    is_survey   = c(0, 0, 1, 1),  # is this a survey (1) or fishery (0)
    qs          = c(1, 1, 6.41, 0.85), # catchability coefficient (q) for surveys
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

# Simulation 1
# Recruitment split equally amongst regions
model_options1 <- model_options
model_options1$recruit_apportionment <- matrix(1/nregions, nrow=nyears+1, ncol=nregions, byrow = TRUE)
model_options1$recruit_apportionment_random <- TRUE
om1 <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = as.matrix(recruitment), 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options1
)

# Simulation 2
# Recruitment split based on proportions from Fenske 2022
model_options2 <- model_options
model_options2$recruit_apportionment <- matrix(rep(c(0.25, 0.43, 0.13, 0.07, 0.14), nyears+1), ncol=nregions, byrow=TRUE)
model_options2$recruit_apportionment_random <- TRUE
om2 <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options2
)

# Simulation 3
# Recruitment split based on the regional proportion of age-2 individuals
model_options3 <- model_options
model_options3$recruit_apportionment <- age2_apportionment
model_options3$recruit_apportionment_random <- TRUE
om3 <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options3
)

# Simulation 4
# Recruitment split based on the regional proportion of spawning biomass
model_options4 <- model_options
model_options4$recruit_apportionment <- ssb_recuit_apportionment
model_options4$recruit_apportionment_random <- TRUE
om4 <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options4
)


# om_sim <- project_multi(
#     init_naa = init_naa, 
#     removals_timeseries = catch_timeseries,
#     recruitment = recruitment, 
#     dem_params = dem_params, 
#     nyears = nyears, 
#     model_options = model_options
# )

# dimnames(om_sim$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("EGOA", "CGOA", "WGOA", "BS", "AI"))
dimnames(om1$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("EGOA", "CGOA", "WGOA", "BS", "AI"))
dimnames(om2$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("EGOA", "CGOA", "WGOA", "BS", "AI"))
dimnames(om3$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("EGOA", "CGOA", "WGOA", "BS", "AI"))
dimnames(om4$naa) <- list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("EGOA", "CGOA", "WGOA", "BS", "AI"))

bind_rows(
    reshape2::melt(om1$naa) %>% mutate(om="om1"),
    reshape2::melt(om2$naa) %>% mutate(om="om2"),
    reshape2::melt(om3$naa) %>% mutate(om="om3"),
    reshape2::melt(om4$naa) %>% mutate(om="om4")
) %>%
    mutate(age_group = case_when(
        age %in% 2:7 ~ "Small",
        age %in% 8:14 ~ "Medium",
        age %in% 15:100 ~ "Large"
    )) %>%
    mutate(
        age_group = factor(age_group)
    ) %>%
    filter(time < 100) %>%
    group_by(time, age_group, region, om) %>%
    summarise(value = mean(sum(value))) %>%
    ggplot(aes(x=time, y=value, fill=age_group)) +
        geom_col(position="fill")+
        scale_fill_manual(values=c("blue", "red", "black"))+
        facet_grid(rows=vars(region), cols=vars(om))+
        theme_bw()+
        coord_cartesian(expand=0)

bind_rows(
    reshape2::melt(apply(om1$naa[1:(nyears),,1,]*dem_params$waa[,,1,], c(1, 2, 3), sum)) %>% mutate(om="om1"),
    reshape2::melt(apply(om2$naa[1:(nyears),,1,]*dem_params$waa[,,1,], c(1, 2, 3), sum)) %>% mutate(om="om2"),
    reshape2::melt(apply(om3$naa[1:(nyears),,1,]*dem_params$waa[,,1,], c(1, 2, 3), sum)) %>% mutate(om="om3"),
    reshape2::melt(apply(om4$naa[1:(nyears),,1,]*dem_params$waa[,,1,], c(1, 2, 3), sum)) %>% mutate(om="om4")
) %>%
    group_by(time, region, om) %>%
    summarise(value = sum(value)) %>%

    ggplot(aes(x=time, y=value, color=om))+
        geom_line()+
        facet_wrap(~region)


m <- bind_rows(
    reshape2::melt(om1$naa) %>% mutate(om="om1"),
    reshape2::melt(om2$naa) %>% mutate(om="om2"),
    reshape2::melt(om3$naa) %>% mutate(om="om3"),
    reshape2::melt(om4$naa) %>% mutate(om="om4")
) %>%
    mutate(age_group = case_when(
        age %in% 2:7 ~ "Small",
        age %in% 8:14 ~ "Medium",
        age %in% 15:100 ~ "Large"
    )) %>%
    filter(time < 100) %>%
    group_by(age_group, region) %>%
    summarise(value = mean(sum(value))) %>%
    pivot_wider(names_from="region", values_from="value") %>%
    ungroup() %>%
    select(EGOA:AI) %>%
    as.matrix

apply(m, 1, \(x) x/sum(x))
