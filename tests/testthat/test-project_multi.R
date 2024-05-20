test_that("pojrect_multi reproduces sablefish timeseries", {
    assessment <- dget("data/test.rdat")

    dem_params <- readRDS("data/sablefish_dem_params.RDS")

    dims <- get_model_dimensions(dem_params$sel)

    init_naa <- array(NA, dim=c(1, dims$nages, dims$nsexes, dims$nregions), dimnames = list(1, 2:31, c("F", "M"), "alaska"))
    init_naa[,,1,] <- assessment$natage.female["1960",]
    init_naa[,,2,] <- assessment$natage.male["1960",]

    recruitment <- assessment$natage.female[,1]*2
    recruitment <- c(recruitment, recruitment[64])

    f_timeseries <- assessment$t.series[,c("F_HAL", "F_TWL")] %>% as.matrix
    f_timeseries <- array(f_timeseries, dim=c(dims$nyears, 1, 1, 1, 2), 
                    dimnames = list("time"=1:dims$nyears, 
                                    age="all",  
                                    sex="all", 
                                    "region"="alaska", 
                                    "fleet"=c("Fixed", "Trawl")))

    model_options <- list(
        removals_input = "F",
        simulate_observations = FALSE
    )

    om_sim <- project_multi(init_naa, f_timeseries, recruitment, dem_params, nyears, model_options)

    om_ssb <- as.vector(apply(om_sim$naa[1:64,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum))
    true_ssb <- assessment$t.series[,"spbiom"]

    expect_equal(om_ssb, true_ssb, tolerance=0.001)

})