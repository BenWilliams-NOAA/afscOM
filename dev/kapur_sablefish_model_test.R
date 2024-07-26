rm(list=ls())
library(devtools)
library(tidyverse)
library(patchwork)

devtools::load_all()

load("data/kapur2024_sablefish_data.rda")

nyears <- length(1960:2019)+1
nages <- length(1:71)
nsexes <- 2
nregions <- 6
nfleets <- 7 

dimension_names <- list(
    "time"=1960:2019,
    "age" = 1:71,
    "sex" = c("F", "M"),
    "region" = c("R1", "R2", "R3", "R4", "R5", "R6"),
    "fleet" = c("F1", "F2", "F3", "F4", "F5", "F6", "F7")
)

region_stock_matrix <- matrix(
    c(
        1,0,0,0,
        0,1,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,1,0,
        0,0,0,1
    ),
    nrow=nregions,
    ncol=4,
    byrow=TRUE
)
dimnames(region_stock_matrix) <- list("region"=c("R1", "R2", "R3", "R4", "R5", "R6"), "stock"=c("S1", "S2", "S3", "S4"))
# t(region_stock_matrix)

region_fleet_matrix <- matrix(
    c(
        0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 1,
        0, 0, 1, 1, 1, 0, 0,
        0, 0, 1, 1, 1, 0, 0,
        1, 1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0    
    ), 
    nrow=nregions, ncol=nfleets,
    byrow=TRUE
)
dimnames(region_fleet_matrix) <- dimension_names[c("region", "fleet")]

###############
## Demographic Parameters
###############
mort_by_stock <- kapur2024_sablefish_data$mortality
mort_by_region <- c(t(mort_by_stock) %*% t(region_stock_matrix)) # convert stock mortality to regional mortalities
mort_matrix <- generate_param_matrix(mort_by_region, dimension_names = dimension_names, by=c("region"))

selex <- aperm(kapur2024_sablefish_data$selectivity[,1:nages,,,drop=FALSE], c(1, 2, 4, 3))
dimnames(selex) <- dimension_names[c("time", "age", "sex", "fleet")]
selex_matrix <- generate_param_matrix(selex, dimension_names = dimension_names, by=c("time", "age", "sex", "fleet"))

# Turn off selectivity for these fleets in these regions
# to indicate which fleets are active in which regions
selex_matrix[,,,1:2, 1:5] <- 0
selex_matrix[,,,3:4, c(1:2, 6:7)] <- 0
selex_matrix[,,,5:6, 3:7] <- 0

movement <- kapur2024_sablefish_data$movement

maturity_by_stock <- kapur2024_sablefish_data$maturity
maturity_by_region <- maturity_by_stock %*% t(region_stock_matrix) # convert stock maturity to regional maturities
maturity_matrix <- generate_param_matrix(maturity_by_region, dimension_names = dimension_names, by=c("age", "region"))

prop_male <- 0.50
sexrat_matrix <- generate_param_matrix(prop_male, dimension_names = dimension_names)

# Kapur et al. 2024 assumes 100% discard mortality for all fleets
dmr_by_region <- rep(1, 6)
names(dmr_by_region) <- colnames(movement)
dmr_matrix <- generate_param_matrix(dmr_by_region, dimension_names = dimension_names, by=c("region"), include_fleet_dim = TRUE)

# Kapur et al. 2024 assumes full retetion for all fleets
retention_matrix <- generate_param_matrix(1, dimension_names = dimension_names, by=c("time", "age", "sex", "fleet"))

waa <- kapur2024_sablefish_data$waa %>% as_tibble() %>% 
    select(sex, age, subarea, weight=mean_weight_at_age_kg)
waa_matrix <- array(NA, dim=c(nyears, nages, nsexes, nregions))

waa_matrix[,,1,1] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R1") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,1] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R1") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,1,2] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R2") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,2] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R2") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,1,3] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R3") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,3] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R3") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,1,4] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R4") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,4] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R4") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,1,5] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R5") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,5] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R5") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,1,6] <- matrix(rep(waa %>% filter(sex == 1, subarea=="R6") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 
waa_matrix[,,2,6] <- matrix(rep(waa %>% filter(sex == 2, subarea=="R6") %>% pull(weight), nyears), nrow=nyears, byrow=TRUE) 


dem_params <- list(
    waa=waa_matrix,
    mat=maturity_matrix,
    mort=mort_matrix,
    sexrat=sexrat_matrix,
    sel=selex_matrix,
    ret=retention_matrix,
    dmr=dmr_matrix,
    surv_sel=NA, # not needed as we aren't going to simulate observations
    movement=movement
)


#######
naa <- aperm(kapur2024_sablefish_data$naa, perm=c(1, 2, 4, 3))
init_naa <- naa[1,,,,drop=FALSE]

nyais <- naa

# Get recruitment from NAA and figure out what the annual
# recruitment to each region is
recruitment <- apply(naa[,1,,], 1, sum)
recruits_by_region <- apply(naa[,1,,], c(1, 3), sum)
recruit_apportionment <- t(apply(recruits_by_region, 1, \(x) x/sum(x)))

# Get catch timeseries by region from catch-at-age data
catch_matrix <- kapur2024_sablefish_data$catch %>% as_tibble() %>%
    group_by(year, fleet) %>%
    summarise(catch = sum(catch_t)) %>%
    pull(catch) %>%
    matrix(., ncol=nfleets, byrow=TRUE)
colnames(catch_matrix) <- c("F1", "F2", "F3", "F4", "F5", "F6", "F7")
catch_timeseries <- apply(catch_matrix, 1, sum)

# Calculate region and region-fleet apportionments scheme as being
# proportional to fleet selected biomass in each region 
selected_biomass <- array(NA, c(60, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
for(i in 1:nfleets){
    selected_biomass[,,,,i] <- nyais[1:60,,,]*waa_matrix[1:60,,,]*selex_matrix[,,,,i]
}

selected_biomass_f <- apply(selected_biomass, c(1, 5), sum)
selected_biomass_rf <- apply(selected_biomass, c(1, 4, 5), sum)

app_rf <- array(NA, dim=c(nyears, nregions, nfleets))
for(y in 1:(nyears-1)){
    for(f in 1:nfleets){
        rat <- selected_biomass_rf[y,,f]/selected_biomass_f[y,f]
        app_rf[y,,f] <- catch_matrix[y,f]*rat
    }
}
app_rf[is.na(app_rf) | is.nan(app_rf) | is.infinite(app_rf)] <- 0

app_rf_perm <- aperm(app_rf, perm=c(1, 3, 2))
for(r in 1:nregions){
    app_rf_perm[,,r] <- t(apply(app_rf_perm[,,r], 1, \(x) x/sum(x)))
}

region_apportionment <- t(apply(apply(app_rf, c(1, 2), sum), 1, \(x) x/sum(x)))
fleet_apportionment <- app_rf_perm
fleet_apportionment[is.nan(fleet_apportionment)] <- 0

regional_catch_timeseries <- sweep(region_apportionment[1:60,], 1, matrix(catch_timeseries, ncol=1), FUN="*")
region_fleet_catch_timeseries <- sweep(fleet_apportionment[1:60,,], c(1, 3), regional_catch_timeseries, FUN="*")

# Set all model options
model_options <- list()
model_options$removals_input = "catch"
model_options$region_apportionment = region_apportionment
model_options$fleet_apportionment = fleet_apportionment
model_options$recruit_apportionment = recruit_apportionment
model_options$recruit_apportionment_random = FALSE
model_options$do_recruits_move = FALSE
model_options$simulate_observations = FALSE

# Run OM simulation
om_sim <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = as.matrix(catch_timeseries),
    recruitment = as.matrix(recruitment), 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

# Plot SSBs from Kapur et al. 2024 and afscOM
ssb_yi <- kapur2024_sablefish_data$ssb
dimnames(ssb_yi) <- dimension_names[c("time", "region")]

kapur2024_bio_plot <- reshape2::melt(ssb_yi) %>%
    ggplot(aes(x=time, y=value, color=region))+
        geom_line(linewidth=1.2)+
        labs(y="SSB", title="Kapur et al. 2024 Regional SSB")+
        coord_cartesian(ylim=c(0, 3e8))+
        theme_bw()

ssbaa <- om_sim$naa[1:(nrow(om_sim$naa)-1),,1,,drop=FALSE]*dem_params$waa[1:(nrow(om_sim$naa)-1),,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE]
afscom_bio_plot <- reshape2::melt(apply(ssbaa, c(1, 4), sum)) %>%
    ggplot(aes(x=time, y=value, color=region))+
        geom_line(linewidth=1.2)+
        labs(y="SSB", title="afscOM Regional SSB")+
        coord_cartesian(ylim=c(0, 3e8))+
        theme_bw()

library(patchwork)

kapur2024_bio_plot + afscom_bio_plot + plot_layout(guides="collect")
