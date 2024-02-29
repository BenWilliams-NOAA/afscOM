test_that("1yr projection w/o movement", {
  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  dem_params <- subset_dem_params(dem_params, r=1, d=4, drop=FALSE) # Make it 1 region

  y <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)

  assessment <- dget("data/test.rdat")
  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]*1e6
  naa[,,2,] <- assessment$natage.male["1960",]*1e6
  
  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1] *1e6 * dem_params$sexrat[,1,,]

  tac <- 3114310
  suppressWarnings({
    out_vars <- project(TAC=tac, dem_params=dem_params, prev_naa=naa, recruitment=rec, fleet.props = c(1, 0), options=model_options)
  })

  total_catch <- apply(out_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-6)

  fleet_catch <- c(apply(out_vars$caa_tmp, c(1, 5), sum))
  expect_equal(fleet_catch, c(tac, 0), tolerance=1e-6)

  ssb <- sum(out_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])/1e6
  expect_equal(ssb, 279.6146, tolerance=1e-4)
})
