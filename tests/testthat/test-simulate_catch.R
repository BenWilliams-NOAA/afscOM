dem_params <- afscOM::sablefish_dem_params
simple_om_spatial <- afscOM::simple_om_spatial
assessment <- afscOM::sablefish_assessment_data

test_that("single year catch simulation for one fleet", {


  model_params <- get_model_dimensions(dem_params$sel)

  fleet_apportionment <- array(
    matrix(c(1, 0), nrow=model_params$nyears, ncol=model_params$nfleets, byrow=TRUE),
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
  naa[,,1,] <- assessment$natage.female["2023",]*1e6
  naa[,,2,] <- assessment$natage.male["2023",]*1e6

  tac <- array(7000, c(64, 1, 1))
  catch <- apportion_catch(
      catch_timeseries = tac,
      apportionment = model_options$fleet_apportionment,
      nyears = model_params$nyears,
      nfleets = model_params$nfleets,
      nregions = model_params$nregions
  )$full_catch

  removals <- subset_matrix(catch, 1, d=1, drop=FALSE)
  suppressWarnings({
    catch_vars <- simulate_catch(
        removals = subset_matrix(removals, r, d=3, drop=TRUE),
        dem_params = dem_params,
        naa = naa,
        options = model_options
    )
  })

  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)
  names(dim(fleet_catch)) <- NULL

  expect_equal(total_catch, 7000, tolerance=1e-2)
  expect_equal(fleet_catch, matrix(c(7000, 0), nrow=1), tolerance=1e-4)
})

test_that("single year catch simulation for two fleets", {

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
  naa[,,1,] <- assessment$natage.female["2023",]*1e6
  naa[,,2,] <- assessment$natage.male["2023",]*1e6

  tac <- array(7000, c(64, 1, 1))
  catch <- apportion_catch(
      catch_timeseries = tac,
      apportionment = model_options$fleet_apportionment,
      nyears = model_params$nyears,
      nfleets = model_params$nfleets,
      nregions = model_params$nregions
  )$full_catch

  removals <- subset_matrix(catch, 1, d=1, drop=FALSE)
  suppressWarnings({
    catch_vars <- simulate_catch(
        removals = subset_matrix(removals, r, d=3, drop=TRUE),
        dem_params = dem_params,
        naa = naa,
        options = model_options
    )
  })

  total_catch <- apply(catch_vars$caa_tmp, 1, sum)
  fleet_catch <- apply(catch_vars$caa_tmp, c(1, 5), sum)
  names(dim(fleet_catch)) <- NULL

  expect_equal(total_catch, 7000, tolerance=1e-2)
  expect_equal(fleet_catch, matrix(c(4900, 2100), nrow=1), tolerance=1e-4)
})


test_that("single year catch simulation with two regions and one fleet", {

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

  caa_tmp = array(NA, dim=c(1, model_params$nages, model_params$nsexes, model_params$nregions, model_params$nfleets))

  tac <- array(100, c(25, 1, 1))
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
      caa_tmp[,,,r,] <- catch_vars$caa_tmp
    }

  })

  total_catch <- apply(caa_tmp, 1, sum)
  region_catch <- apply(caa_tmp, c(1, 3), sum)

  expect_equal(total_catch, 100, tolerance=1e-4)
  expect_equal(apply(caa_tmp, c(1, 4), sum), matrix(c(50, 50), ncol=2), tolerance=1e-4)

})

test_that("Apportion catch with recruitment props", {
    nyears <- 10
    nregions <- 5
    nfleets <- 4
    catch_timeseries <- rep(20, nyears)

    fleet_apportionment <- array(
      matrix(c(0.5, 0.3, 0.1, 0.1), nrow=nyears, ncol=nfleets, byrow=TRUE),
      dim=c(nyears, nfleets, nregions)
    )

    model_options <- list(
      fleet_apportionment = fleet_apportionment,
      removals_input = "catch",
      simulate_observations = FALSE
    )

    c <- apportion_catch(
      catch_timeseries = catch_timeseries,
      apportionment = model_options$fleet_apportionment,
      nyears = nyears,
      nfleets = nfleets,
      nregions = nregions
    )

    true_catch <- sweep(fleet_apportionment, 1, catch_timeseries, FUN="*")
    expect_equal(c$full_catch, true_catch)
})


test_that("Apportion catch with catch matrix", {
    nyears <- 10
    nregions <- 5
    nfleets <- 4
    catch_timeseries <- rep(20, nyears)

    fleet_apportionment <- array(
      matrix(c(0.5, 0.3, 0.1, 0.1)/nregions, nrow=nyears, ncol=nfleets, byrow=TRUE),
      dim=c(nyears, nfleets, nregions)
    )

    true_catch <- sweep(fleet_apportionment, 1, catch_timeseries, FUN="*")

    model_options <- list(
      removals_input = "catch",
      simulate_observations = FALSE
    )

    c <- apportion_catch(
      catch_timeseries = true_catch,
      apportionment = model_options$fleet_apportionment,
      nyears = nyears,
      nfleets = nfleets,
      nregions = nregions
    )

    expect_equal(c$full_catch, true_catch)
    expect_equal(c$region_fleet_props, fleet_apportionment, tolerance=1e-5)
})
