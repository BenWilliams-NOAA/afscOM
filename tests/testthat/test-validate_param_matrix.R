dem_params <- afscOM::sablefish_dem_params


test_that("Missing DMR", {

    model_params <- list(nyears=64, nages=30, nsexes=2, nregions=1, nfleets=2)



    expect_equal(validate_dem_params(dem_params, model_params), dem_params)

    dem_params <- dem_params[names(dem_params) != "dmr"]

    expect_error(validate_dem_params(dem_params, model_params))
})

test_that("Wrong dimensions", {

    model_params <- list(nyears=64, nages=30, nsexes=2, nregions=1, nfleets=2)

    expect_equal(validate_dem_params(dem_params, model_params), dem_params)

    dem_params$mort <- dem_params$mort[1:20,1:15,1,1,drop=FALSE]
    expect_error(validate_dem_params(dem_params, model_params))
})
