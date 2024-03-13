test_that("1yr projection w/o movement", {
  
  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(1, 0),
    fleet_apportionment = c(1, 0),
    removals_input = "catch",
    simulate_observations = FALSE
  )
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  assessment <- dget("data/test.rdat")
  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]
  naa[,,2,] <- assessment$natage.male["1960",]
  
  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1]

  tac <- 3.114310
  suppressWarnings({
    out_vars <- project(removals=tac, dem_params=dem_params, prev_naa=naa, recruitment=rec, fleet.props = c(1, 0), options=model_options)
  })

  total_catch <- apply(out_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-2)

  ssb <- sum(out_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])
  expect_equal(ssb, 278.43, tolerance=1e-4)
})

test_that("1yr projection w F input", {
  
  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(1, 0),
    fleet_apportionment = c(1, 0),
    removals_input = "F",
    simulate_observations = FALSE
  )
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  assessment <- dget("data/test.rdat")
  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]
  naa[,,2,] <- assessment$natage.male["1960",]
  
  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1]

  tac <- 3.701
  f_timeseries <- array(c(0.006515313, 0.000000000), dim=c(1, 1, 1, 1, 2))
  f_timeseries <- subset_matrix(f_timeseries, y, d=1, drop=FALSE)
  suppressWarnings({
    out_vars <- project(removals=f_timeseries, dem_params=dem_params, prev_naa=naa, recruitment=rec, fleet.props = c(1, 0), options=model_options)
  })

  total_catch <- apply(out_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-2)

  ssb <- sum(out_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])
  expect_equal(ssb, 278.14, tolerance=1e-4)
})
