rm(list=ls())
library(tidyverse)
library(patchwork)
devtools::load_all()

# Load OM objects for simple OM
load("data/simple_om.rda")

# initial numbers-at-age
simple_om$init_naa

# recruitment timeseries
simple_om$recruitment

# catch timeseries
simple_om$catch

# get model dimensions
model_dimensions <- get_model_dimensions(simple_om$dem_params$sel)

simple_om$options$region_apportionment <- NULL
simple_om$options$fleet_apportionment <- array(1, dim=c(model_dimensions$nyears, model_dimensions$nfleet, model_dimensions$nregions))
simple_om$options$recruit_apportionment = array(1, dim=c(model_dimensions$nyears+1, 1))
simple_om$options$random_apportion_recruits = FALSE
simple_om$options$do_recruits_move = FALSE
simple_om$options$simulate_observations = FALSE
simple_om$options$recruitment_pars = list(
    h = 0.75,
    R0 = median(simple_om$recruitment),
    S0 = 2000
)

# run the om simulatino
om_sim <- project_multi(
    init_naa = simple_om$init_naa, 
    removals_timeseries = simple_om$catch,
    recruitment = beverton_holt, 
    dem_params = simple_om$dem_params, 
    nyears = model_dimensions$nyears, 
    model_options = simple_om$options
)

dimnames(om_sim$naa) <- list("time"=1:(model_dimensions$nyears+1), "age"=1:10, "sex"=c("F", "M"), "region"="Region 1")
dimnames(om_sim$caa) <- list("time"=1:(model_dimensions$nyears), "age"=1:10, "sex"=c("F", "M"), "region"="Region 1", "fleet"=("Fleet 1"))


om_sim$recruits

# derived quantities and plots
ssb <- compute_ssb(om_sim$naa, simple_om$dem_params)
bio <- compute_bio(om_sim$naa, simple_om$dem_params)
catch <- compute_total_catch(om_sim$caa)
f <- compute_total_f(om_sim$faa)

ssb_plot <- plot_ssb(ssb)
bio_plot <- plot_bio(bio)
catch_plot <- plot_catch(catch)
f_plot <- plot_f(f)

((ssb_plot + bio_plot) / (catch_plot + f_plot)) + plot_layout(guides="collect")

naa_plot <- plot_atage(om_sim$naa)
caa_plot <- plot_atage(om_sim$caa)

naa_plot
caa_plot

