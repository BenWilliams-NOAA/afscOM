dem_params <- afscOM::sablefish_dem_params
assessment <- afscOM::sablefish_assessment_data

test_that("1yr projection w/o movement", {

  model_params <- get_model_dimensions(dem_params$sel)

  fleet_apportionment <- array(
    matrix(c(1, 0), nrow=model_params$nyears, ncol=model_params$nfleets, byrow=TRUE),
    dim=c(model_params$nyears, model_params$nfleets, model_params$nregions)
  )

  model_options <- list(
    # regional_apportionment = matrix(1, nrow = model_params$nyears, ncol=model_params$nregions),
    fleet_apportionment = fleet_apportionment,
    removals_input = "catch",
    simulate_observations = FALSE
  )

  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]
  naa[,,2,] <- assessment$natage.male["1960",]

  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1]

  tac <- 3.114310
  removals <- apportion_catch(
      catch_timeseries = tac,
      apportionment = model_options$fleet_apportionment,
      nyears = model_params$nyears,
      nfleets = model_params$nfleets,
      nregions = model_params$nregions
  )$full_catch

  removals_input <- subset_matrix(removals, 1, d=1, drop=FALSE)

  suppressWarnings({
    out_vars <- project_single(
        removals=removals_input,
        dem_params=dem_params,
        prev_naa=naa,
        recruitment=rec,
        options=model_options
    )
  })

  total_catch <- apply(out_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-2)

  ssb <- sum(out_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])
  expect_equal(ssb, 281, tolerance=0.001)
})

test_that("1yr projection w F input", {

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

  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]
  naa[,,2,] <- assessment$natage.male["1960",]

  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1]

  tac <- 3.114
  f_timeseries <- array(c(0.006515313, 0.000000000), dim=c(1, 2, 1))
  f_timeseries <- subset_matrix(f_timeseries, y, d=1, drop=FALSE)
  suppressWarnings({
    out_vars <- project_single(removals=f_timeseries, dem_params=dem_params, prev_naa=naa, recruitment=rec, options=model_options)
  })

  total_catch <- apply(out_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-2)

  ssb <- sum(out_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])
  expect_equal(ssb, 281, tolerance=1e-3)
})
