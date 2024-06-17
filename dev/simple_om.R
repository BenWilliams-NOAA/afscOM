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

# run the om simulatino
om_sim <- project_multi(
    init_naa = simple_om$init_naa, 
    removals_timeseries = simple_om$catch,
    recruitment = simple_om$recruitment, 
    dem_params = simple_om$dem_params, 
    nyears = model_dimensions$nyears, 
    model_options = simple_om$options
)

dimnames(om_sim$naa) <- list("time"=1:(model_dimensions$nyears+1), "age"=1:10, "sex"=c("F", "M"), "region"="Region 1")
dimnames(om_sim$caa) <- list("time"=1:(model_dimensions$nyears), "age"=1:10, "sex"=c("F", "M"), "region"="Region 1", "fleet"=("Fleet 1"))


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

