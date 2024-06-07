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

# derived quantities and plots
ssb <- apply(om_sim$naa[1:25,,1,]*simple_om$dem_params$waa[,,1,]*simple_om$dem_params$mat[,,1,], 1, sum)
bio <- apply(om_sim$naa[1:25,,,]*simple_om$dem_params$waa[,,,], 1, sum)
catch <- apply(om_sim$caa, 1, sum)
f <- apply(apply(om_sim$faa, c(1, 5), \(x) max(x)), 1, sum)

ssb_comp <- data.frame(
    year=1:25,
    om_ssb=ssb,
    om_catch=catch,
    om_bio = bio,
    om_f = f
)

p1 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=om_ssb, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 7000), breaks=seq(0, 7000, 1000))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="SSB", x="Year", title="Spawning Biomass Comparison")+
    theme_bw()

p2 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=om_catch, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 25))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="Catch", x="Year", title="Total Catch Comparison")+
    theme_bw()

p3 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=om_bio, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 15000), breaks=seq(0, 15000, 2000))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="Biomass", x="Year", title="Total Biomass Comparison")+
    theme_bw()

p4 <- ggplot(ssb_comp, aes(x=year))+
    geom_line(aes(y=om_f, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 0.04), breaks=seq(0, 0.04, 0.01))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="F", x="Year", title="Fishing Mortality Comparison")+
    theme_bw()

((p1+p3) / (p2+p4)) + plot_layout(guides="collect")
