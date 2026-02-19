assessment <- afscOM::sablefish_assessment_data
dem_params <- afscOM::sablefish_dem_params
simple_om_spatial <- afscOM::simple_om_spatial

test_that("single year age-structure population", {


  model_params <- get_model_dimensions(dem_params$sel)

  fleet_apportionment <- array(
    matrix(c(0.7, 0.3), nrow=model_params$nyears, ncol=model_params$nfleets, byrow=TRUE),
    dim=c(model_params$nyears, model_params$nfleets, model_params$nregions)
  )

  model_options <- list(
    regional_apportionment = 1,
    fleet_apportionment = fleet_apportionment,
    removals_input = "catch"
  )

  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)
  dem_params <- subset_dem_params(dem_params, r, d=4, drop=FALSE)

  naa <- array(NA, dim=c(1, model_params$nages, model_params$nsexes, 1))
  naa[,,1,] <- assessment$natage.female["1960",]*1e6
  naa[,,2,] <- assessment$natage.male["1960",]*1e6

  rec <- array(NA, dim=c(1, 1, model_params$nsexes, 1))
  rec[1,1,,1] <- 2*assessment$natage.female["1960",1] *1e6 * dem_params$sexrat[,1,,]

  tac <- 3114310
  suppressWarnings({

    catch <- apportion_catch(
        catch_timeseries = tac,
        apportionment = model_options$fleet_apportionment,
        nyears = model_params$nyears,
        nfleets = model_params$nfleets,
        nregions = model_params$nregions
    )$full_catch

    removals <- subset_matrix(catch, 1, d=1, drop=FALSE)

    catch_vars <- simulate_catch(subset_matrix(removals, r, d=3, drop=TRUE), dem_params=dem_params, naa=naa, options=model_options)
    pop_vars <- simulate_population(prev_naa=naa, faa=catch_vars$faa_tmp, recruitment=rec, dem_params=dem_params, options=options)
  })

  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-6)

  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)
  names(dim(fleet_catch)) <- NULL
  expect_equal(fleet_catch, matrix(c(2180017, 934293), nrow=1), tolerance=1e-5)

  ssb <- sum(pop_vars$naa[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,])/1e6
  expect_equal(ssb, 281, tolerance=1e-2)
})


test_that("single year population with two regions", {

  dem_params <- simple_om_spatial$dem_params
  model_params <- get_model_dimensions(dem_params$sel)
  model_options <- simple_om_spatial$model_options

  y <- 1
  dem_params <- subset_dem_params(dem_params, y, d=1, drop=FALSE)

  fleet_apportionment <- aperm(array(
    matrix(c(0.5, 0.5), nrow=model_params$nyears, ncol=model_params$nregions, byrow=TRUE),
    dim=c(model_params$nyears, model_params$nregions, model_params$nfleets)
  ), c(1, 3, 2))

  model_options$fleet_apportionment <- fleet_apportionment


  caa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))
  naa_tmp         = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions))

  tac <- 100
  suppressWarnings({

    catch <- apportion_catch(
        catch_timeseries = tac,
        apportionment = model_options$fleet_apportionment,
        nyears = model_params$nyears,
        nfleets = model_params$nfleets,
        nregions = model_params$nregions
    )$full_catch

    removals <- subset_matrix(catch, 1, d=1, drop=FALSE)

    for(r in 1:model_params$nregions){
      remove <- subset_matrix(removals, r, d=3, drop=TRUE)
      dp.r <- subset_dem_params(dem_params, r, d=4, drop=FALSE)
      naa.r <- subset_matrix(simple_om_spatial$init_naa, r=r, d=4, drop=FALSE)
      catch_vars <- simulate_catch(
        removals=remove,
        dem_params=dp.r,
        naa=naa.r,
        options=model_options
      )
      pop_vars <- simulate_population(prev_naa=naa.r, faa=catch_vars$faa_tmp, recruitment=60, dem_params=dp.r, options=options)
      caa_tmp[,,,r,] <- catch_vars$caa_tmp
      naa_tmp[,,,r] <- pop_vars$naa
    }
  })

  total_catch <- apply(caa_tmp, 1, sum)
  expect_equal(total_catch, tac, tolerance=1e-6)

  region_catch <- apply(caa_tmp, c(1, 4), sum)
  expect_equal(region_catch, matrix(c(50, 50), nrow=1), tolerance=1e-5)

  ssb <- apply(naa_tmp[,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 2, sum)
  expect_equal(as.vector(ssb), rep(2079.496, 2), tolerance=1e-2)

})
