test_that("single year age-structure population", {

  load(file.path(here::here(), "data/sablefish_dem_params.rda"))
  dem_params <- sablefish_dem_params
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(0.70, 0.30),
    fleet_apportionment = c(0.70, 0.30),
    removals_input = "catch"
  )
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
  assessment <- sablefish_assessment_data
  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]*1e6
  naa[,,2,] <- assessment$natage.male["1960",]*1e6
  
  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1] *1e6 * dem_params$sexrat[,1,,]

  tac <- 3114310
  suppressWarnings({
    catch_vars <- simulate_catch(tac, fleet.props=c(0.70, 0.30), dem_params=dem_params, naa=naa, options=model_options)
    pop_vars <- simulate_population(prev_naa=naa, faa=catch_vars$faa_tmp, recruitment=rec, dem_params=dem_params, options=options)
  })

  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-6)

  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)
  expect_equal(fleet_catch, matrix(c(2180017, 934293), nrow=1), tolerance=1e-5)

  ssb <- sum(pop_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])/1e6
  expect_equal(ssb, 281, tolerance=1e-2)
})
