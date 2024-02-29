test_that("single year age-structure population", {

  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(0.70, 0.30),
    fleet_apportionment = c(0.70, 0.30)
  )
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  assessment <- dget("data/test.rdat")
  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]*1e6
  naa[,,2,] <- assessment$natage.male["1960",]*1e6
  
  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1] *1e6 * dem_params$sexrat[,1,,]

  tac <- 3114310
  suppressWarnings({
    catch_vars <- simulate_catch(tac, fleet.props=c(1.0, 0.0), dem_params=dem_params, naa=naa, options=model_options)
    pop_vars <- simulate_population(prev.naa=naa, faa=catch_vars$faa, recruitment=rec, dem_params=dem_params, options=options)
  })

  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-6)

  ssb <- sum(pop_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])/1e6
  expect_equal(ssb, 279.8009, tolerance=1e-4)
})
