test_that("single year catch simulation", {

  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(0.70, 0.30),
    fleet_apportionment = c(0.70, 0.30)
  )
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  assessment <- dget("data/test.rdat")
  naa <- matrix(NA, nrow=model_params$nages, ncol=2)
  naa[,1] <- assessment$natage.female["2023",]*1e6
  naa[,2] <- assessment$natage.male["2023",]*1e6

  tac <- 7000
  suppressWarnings({
    catch_vars <- simulate_catch(tac, fleet.props=c(0.70, 0.30), dem_params=dem_params, naa=naa, options=model_options)
  })
  
  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)

  expect_equal(total_catch, 7000)
  expect_equal(fleet_catch, matrix(c(4900, 2100), nrow=1), tolerance=1e-6)
})
