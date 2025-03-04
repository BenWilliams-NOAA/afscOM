test_that("project_multi reproduces sablefish timeseries", {
    
    load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
    assessment <- sablefish_assessment_data

    load(file.path(here::here(), "data/sablefish_dem_params.rda"))
    dem_params <- sablefish_dem_params

    dims <- get_model_dimensions(dem_params$sel)

    init_naa <- array(NA, dim=c(1, dims$nages, dims$nsexes, dims$nregions), dimnames = list(1, 2:31, c("F", "M"), "alaska"))
    init_naa[,,1,] <- assessment$natage.female["1960",]
    init_naa[,,2,] <- assessment$natage.male["1960",]

    recruitment <- assessment$natage.female[,1]*2
    recruitment <- as.matrix(c(recruitment, recruitment[64]))

    f_timeseries <- assessment$t.series[,c("F_HAL", "F_TWL")] %>% as.matrix
    f_timeseries <- array(f_timeseries, dim=c(dims$nyears, 1, 1, 1, 2), 
                    dimnames = list("time"=1:dims$nyears, 
                                    age="all",  
                                    sex="all", 
                                    "region"="alaska", 
                                    "fleet"=c("Fixed", "Trawl")))

    model_options <- list(
        removals_input = "F",
        simulate_observations = FALSE,
        recruit_apportionment = NULL,
        recruit_apportionment_random = FALSE
    )

    om_sim <- project_multi(
        init_naa, 
        f_timeseries, 
        recruitment, 
        dem_params, 
        nyears, 
        model_options
    )

    om_ssb <- as.vector(apply(om_sim$naa[1:64,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum))
    true_ssb <- assessment$t.series[,"spbiom"]

    expect_equal(om_ssb, true_ssb, tolerance=0.001)

})

test_that("project_multi reproduces sablefish timeseries with catch timeseries", {
    
    load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
    assessment <- sablefish_assessment_data

    load(file.path(here::here(), "data/sablefish_dem_params.rda"))
    dem_params <- sablefish_dem_params

    dims <- get_model_dimensions(dem_params$sel)

    init_naa <- array(NA, dim=c(1, dims$nages, dims$nsexes, dims$nregions), dimnames = list(1, 2:31, c("F", "M"), "alaska"))
    init_naa[,,1,] <- assessment$natage.female["1960",]
    init_naa[,,2,] <- assessment$natage.male["1960",]

    recruitment <- assessment$natage.female[,1]*2
    recruitment <- as.matrix(c(recruitment, recruitment[64]))

    TACs <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
    catch_timeseries <- assessment$t.series[,c("Catch_HAL", "Catch_TWL")] %>% as.matrix
    catch_timeseries <- array(catch_timeseries, dim=c(nyears, nfleets, nregions),
                    dimnames = list("time"=1:nyears,
                                    "fleet"=c("Fixed", "Trawl"),
                                    "region"="alaska")
    )

    model_options <- list(
        removals_input = "catch",
        simulate_observations = FALSE,
        recruit_apportionment = NULL,
        recruit_apportionment_random = FALSE
    )
    suppressWarnings({
        om_sim <- project_multi(
            init_naa, 
            catch_timeseries, 
            recruitment, 
            dem_params, 
            nyears, 
            model_options
        )
    })
    

    om_ssb <- as.vector(apply(om_sim$naa[1:64,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum))
    true_ssb <- assessment$t.series[,"spbiom"]

    expect_equal(om_ssb, true_ssb, tolerance=0.01)

})
