test_that("single year catch simulation for one fleet", {

  load(file.path(here::here(), "data/sablefish_dem_params.rda"))
  dem_params <- sablefish_dem_params
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = 1,
    fleet_apportionment = c(1.00, 0.00),
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
  naa[,,1,] <- assessment$natage.female["2023",]*1e6
  naa[,,2,] <- assessment$natage.male["2023",]*1e6

  tac <- 7000
  suppressWarnings({
    catch_vars <- simulate_catch(tac, fleet_props=matrix(c(1.00, 0.00), ncol=model_params$nfleets), dem_params=dem_params, naa=naa, options=model_options)
  })
  
  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)

  expect_equal(total_catch, 7000)
  expect_equal(fleet_catch, matrix(c(7000, 0), nrow=1), tolerance=1e-4)
})

test_that("single year catch simulation for two fleets", {

  load(file.path(here::here(), "data/sablefish_dem_params.rda"))
  dem_params <- sablefish_dem_params
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- list(
    regional_apportionment = c(1, 0),
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
  naa[,,1,] <- assessment$natage.female["2023",]*1e6
  naa[,,2,] <- assessment$natage.male["2023",]*1e6

  tac <- 7000
  suppressWarnings({
    catch_vars <- simulate_catch(tac, fleet_props=matrix(c(0.70, 0.30), ncol=model_params$nfleets), dem_params=dem_params, naa=naa, options=model_options)
  })
  
  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)

  expect_equal(total_catch, 7000, tolerance=1e-4)
  expect_equal(fleet_catch, matrix(c(4900, 2100), nrow=1), tolerance=1e-4)
})


test_that("single year catch simulation with two regions and one fleet", {
  
  load(file.path(here::here(), "data", "simple_om_spatial.rda"))
  dem_params <- simple_om_spatial$dem_params
  model_params <- get_model_dimensions(dem_params$sel)

  y <- 1
  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  region_props <- model_options$region_apportionment[y,,drop=FALSE]
  fleet_props <- model_options$fleet_apportionment[y,,drop=FALSE]

  caa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))

  tac <- 100
  suppressWarnings({
    for(r in 1:model_params$nregions){
      remove <- tac*region_props[1,r]
      dp.r <- subset_dem_params(dem_params, r, d=4, drop=FALSE)
      naa.r <- subset_matrix(simple_om_spatial$init_naa, r=r, d=4, drop=FALSE)
      catch_vars <- simulate_catch(
        removals=remove, 
        fleet_props=fleet_props, 
        dem_params=dp.r, 
        naa=naa.r, 
        options=model_options
      )
      caa_tmp[,,,r,] <- catch_vars$caa_tmp
    }

    expect_equal(sum(caa_tmp), tac, tolerance=1e-4)
    expect_equal(apply(caa_tmp, c(1, 4), sum), matrix(c(50, 50), ncol=2), tolerance=1e-4)
  })



})