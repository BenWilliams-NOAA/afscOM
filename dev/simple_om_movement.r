nyears <- 25
nages <- 10
nsexes <- 2
nregions <- 2
nfleets <- 1

dimension_names = list("time"=1:nyears, "age"=1:nages, "sex"=c("F", "M"), "region"=c("Region 1", "Region 2"), "fleet"=c("Fleet 1"))

M <- 0.20
mort <- generate_param_matrix(M, dimension_names = dimension_names)

prop_males <- 0.5
sexrat <- generate_param_matrix(prop_males, dimension_names = dimension_names)

retention <- 1.0
ret <- generate_param_matrix(retention, dimension_names = dimension_names, include_fleet_dim = TRUE)

discard <- 0.0
dmr <- generate_param_matrix(discard, dimension_names = dimension_names, include_fleet_dim = TRUE)

maturity <- c(0, 0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
mat <- generate_param_matrix(maturity, dimension_names = dimension_names, by="age")

weight_mat <- matrix(NA, nrow=nages, ncol=nsexes, dimnames=dimension_names[c("age", "sex")])
weight_mat[,1] <- 1:10
weight_mat[,2] <- 0.80*weight_mat[,1]
waa <- generate_param_matrix(weight_mat, dimension_names = dimension_names, by=c("age", "sex"))

selex_mat <- array(NA, dim=c(nages, nsexes, nregions, nfleets), dimnames=dimension_names[c("age", "sex", "region", "fleet")])
selex_mat[,1,1,1] <- c(0, 0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
selex_mat[,2,1,1] <- c(0, 0.20, 0.40, 0.65, 0.90, 1.0, 1.0, 1.0, 1.0, 1.0)
selex_mat[,1,2,1] <- c(0, 0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
selex_mat[,2,2,1] <- c(0, 0.20, 0.40, 0.65, 0.90, 1.0, 1.0, 1.0, 1.0, 1.0)
sel <- generate_param_matrix(selex_mat, dimension_names = dimension_names, by=c("age", "sex", "region", "fleet"), include_fleet_dim = TRUE)

survey_selex_mat <- array(NA, dim=c(nages, nsexes, nregions, nfleets), dimnames=dimension_names[c("age", "sex", "region", "fleet")])
survey_selex_mat[,1,1,1] <- c(0, 0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
survey_selex_mat[,2,1,1] <- c(0, 0.20, 0.40, 0.65, 0.90, 1.0, 1.0, 1.0, 1.0, 1.0)
survey_selex_mat[,1,2,1] <- c(0, 0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
survey_selex_mat[,2,2,1] <- c(0, 0.20, 0.40, 0.65, 0.90, 1.0, 1.0, 1.0, 1.0, 1.0)
survey_sel <- generate_param_matrix(survey_selex_mat, dimension_names = dimension_names, by=c("age", "sex", "region", "fleet"), include_fleet_dim = TRUE)

young_move <- matrix(c(0.2, 0.8, 0.05, 0.95), ncol=nregions, byrow=TRUE)
old_move <- matrix(c(0.99, 0.01, 0.15, 0.85), ncol=nregions, byrow=TRUE)



movement_matrix <- array(NA, dim=c(nregions, nregions, nages))
movement_matrix[,,1:5] <- young_move
movement_matrix[,,6:10] <- old_move

dem_params <- list(
    waa=waa,
    mat=mat,
    mort=mort,
    sexrat=sexrat,
    sel=sel,
    ret=ret,
    dmr=dmr,
    surv_sel=survey_sel,
    movement = movement_matrix
)

init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames = list(1, 1:10, c("F", "M"), c("Region 1", "Region 2")))
init_naa[,,1,] <- c(100, 90, 80, 70, 60, 50, 40, 30, 20, 30)
init_naa[,,2,] <- c(100, 90, 80, 70, 60, 50, 40, 30, 20, 30)

recruitment <- 2*c(30, 30, 30, 30, 30, 30, 5000, 60, 30, 30, 30, 30, 30, 50, 70, 30, 30, 30, 5000, 30, 30, 30, 30, 30, 30)

catch_timeseries <- as.matrix(rep(50, nyears))

model_options <- list()
model_options$removals_input = "catch"
model_options$region_apportionment = matrix(rep(c(0.50, 0.50), nyears), ncol=nregions)
model_options$fleet_apportionment = matrix(1, nrow=nyears, ncol=nfleets)
model_options$simulate_observations = FALSE

simple_om_spatial <- list(
    dem_params = dem_params,
    init_naa = init_naa,
    recruitment = recruitment,
    catch = catch_timeseries,
    model_options=model_options
)
usethis::use_data(simple_om_spatial, overwrite = TRUE)

om_sim <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = catch_timeseries,
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

devtools::load_all()


ssb <- apply(om_sim$naa[1:25,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], c(1, 3), sum)
catch <- apply(om_sim$caa, c(1, 4), sum)

plot(1:25, ssb[,1], type="l", ylim=c(0, 12000))
lines(1:25, ssb[,2], col="red")






f_timeseries <- om_sim$f

model_options <- list()
model_options$removals_input = "F"
# model_options$region_apportionment = matrix(rep(c(0.50, 0.50), nyears), ncol=nregions)
# model_options$fleet_apportionment = matrix(1, nrow=nyears, ncol=nfleets)
model_options$simulate_observations = FALSE

om_sim2 <- project_multi(
    init_naa = init_naa, 
    removals_timeseries = f_timeseries,
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

om_sim$naa
om_sim2$naa

om_sim$naa == om_sim2$naa

