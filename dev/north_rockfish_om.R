#' Northern Rockfish OM
#' Based on November 2024 Northern Rockfish
#' stock assessment model (Williams et al. 2024) 


library(devtools)
devtools::load_all()

data <- readRDS("~/Downloads/dat.RDS")
assessment <- readRDS("~/Downloads/m24.RDS")

nyears <- length(data$years)
nages  <- 50
nlengths <- length(data$length_bins)
nsexes <- 1
nregions <- 1
nfleets <- 1
nsurveys <- 1

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 1:50,
    "sex"  = c("F"),
    "region" = "GOA",
    "fleet" = "F1"
)

model_params <- set_model_params(nyears, nages, nlengths, nsexes, nregions, nfleets)
model_options <- setup_model_options(model_params)

#' 1. Setup demographic rates
M <- 0.059
mort_matrix <- generate_param_matrix(M, dimension_names = dimension_names)

waa <- data$waa
waa_matrix <- generate_param_matrix(waa, dimension_names = dimension_names, by=c("age"))

mat <- data$maa
mat_matrix <- generate_param_matrix(mat, dimension_names = dimension_names, by="age")

sexrat <- 1
sexrat_matrix <- generate_param_matrix(sexrat, dimension_names = dimension_names)

ret <- 1
ret_matrix <- generate_param_matrix(ret, dimension_names = dimension_names, include_fleet_dim = TRUE)

dmr <- 1
dmr_matrix <- generate_param_matrix(dmr, dimension_names = dimension_names, include_fleet_dim = TRUE)

delta <- assessment$deltaC
a50 <- assessment$a50C
sel <- 1/(1+exp(-delta*((1:50)-a50)))
sel_matrix <- generate_param_matrix(sel, dimension_names=dimension_names, by="age", include_fleet_dim = TRUE)

dem_params <- list(
    waa=waa_matrix,
    mat=mat_matrix,
    mort=mort_matrix,
    sexrat=sexrat_matrix,
    sel=sel_matrix,
    ret=ret_matrix,
    dmr=dmr_matrix
)

#' 2. Catch timeseries
removals <- data$catch_obs
removals_matrix <- array(removals, dim=c(nyears, 1, nsexes, nregions, nfleets))

f_matrix <- array(assessment$Ft, dim=c(nyears, 1, nsexes, nregions, nfleets))

#' 3. Initial NAA
init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames = list(1, 1:nages, c("F"), "alaska"))
init_naa[,,1,] <- assessment$Nat[,1]

#' 4. Recruitment
recs <- c(assessment$Nat[1,], mean(assessment$Nat[1,]))

#' 5. Model Options
model_options <- setup_model_options(model_dimensions = model_params)
model_options$random_apportion_recruits <- FALSE
model_options$simulate_observations <- FALSE
model_options$removals_input <- "F"

#' 6. Run OM Simulation
om_sim <- project(
    init_naa = init_naa, 
    removals_timeseries = f_matrix, 
    recruitment = recs, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

#' 7. Outputs
ssb <- compute_ssb(om_sim$naa, dem_params)/2
bio <- compute_bio(om_sim$naa, dem_params)
catch <- compute_total_catch(om_sim$caa)
fleet_catch <- compute_fleet_catch(om_sim$caa)
f <- compute_total_f(om_sim$faa)

par(mfrow=c(2, 2))
plot(1:nyears, assessment$spawn_bio, type="l")
lines(1:nyears, ssb, col="red")
plot(1:nyears, assessment$tot_bio, type="l")
lines(1:nyears, bio, col="red")
plot(1:nyears, assessment$catch_pred, type="l")
lines(1:nyears, catch, col="red")
plot(1:nyears, assessment$Ft, type="l")
lines(1:nyears, f, col="red")
