test_that("Recruitment, single region", {
    nyears <- 10
    nregions <- 1
    recruitment <- rep(20, nyears+1)
    model_options <- list()

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    expect_equal(r$rec_props, array(1, dim=c(11, 1)))
    expect_equal(r$full_recruitment, array(recruitment, dim=c(nyears+1, 1)))
})

test_that("Recruitment vector no apportionment", {
    nyears <- 10
    nregions <- 5
    recruitment <- rep(20, nyears+1)
    model_options <- list()

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    expect_equal(r$rec_props, array(0.2, dim=c(11, 5)))
    expect_equal(r$full_recruitment, array(4, dim=c(11, 5)))

})

test_that("Recruitment array no apportionment", {
    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list()

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    expect_equal(r$rec_props, array(0.2, dim=c(11, 5)))
    expect_equal(r$full_recruitment, array(4, dim=c(11, 5)))

})

test_that("Recruitment array, fixed apportionment", {
    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list(
        recruit_apportionment = c(0.5, 0.2, 0.1, 0.1, 0.1)
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    expect_equal(r$rec_props, matrix(c(0.5, 0.2, 0.1, 0.1, 0.1), nrow=nyears+1, ncol=nregions, byrow=TRUE))
    expect_equal(r$full_recruitment, matrix(c(10, 4, 2, 2, 2), nrow=nyears+1, ncol=nregions, byrow=TRUE))

}) 

test_that("Recruitment array, time-varying apportionment", {
    nyears <- 4
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))

    apportionment_matrix <- matrix(
                                c(
                                    0.5, 0.2, 0.1, 0.1, 0.1,
                                    0.1, 0.5, 0.2, 0.1, 0.1,
                                    0.1, 0.1, 0.5, 0.2, 0.1,
                                    0.1, 0.1, 0.1, 0.5, 0.2,
                                    0.2, 0.1, 0.1, 0.1, 0.5
                                ),
                                nrow=nyears+1,
                                ncol=nregions,
                                byrow=TRUE
                            )

    model_options <- list(
        recruit_apportionment = apportionment_matrix
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    true_recruitment <- sweep(apportionment_matrix, 1, recruitment, FUN="*")
    expect_equal(r$rec_props, apportionment_matrix)
    expect_equal(r$full_recruitment, true_recruitment)

}) 

test_that("Recruitment array, apportionment function", {
    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list(
        recruit_apportionment = mean
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    expect_type(r$rec_props, "closure")
    expect_equal(r$full_recruitment, recruitment)

})


test_that("Get annual recuits, fixed apportionment", {

    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list(
        recruit_apportionment = c(0.5, 0.2, 0.1, 0.1, 0.1),
        recruit_apportionment_random = FALSE
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    r_y <- get_annual_recruitment(
        y = 1,
        rec_timeseries = r$full_recruitment,
        apportionment = r$rec_props,
        apportion_random = model_options$recruit_apportionment_random,
        apportionment_pars = model_options$recruit_apportionment_pars
    )

    expect_equal(r_y, array(20*c(0.5, 0.2, 0.1, 0.1, 0.1), dim=c(1, 5)))

})

test_that("Get annual recuits, function apportionment", {

    example_rec_generation_function <- function(input1, input2){
        return(c(0.5, 0.2, 0.1, 0.1, 0.1))
    }

    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list(
        recruit_apportionment = example_rec_generation_function,
        recruit_apportionment_pars = list(
            input1 = 12,
            input2 = 24
        ),
        recruit_apportionment_random = FALSE
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    r_y <- get_annual_recruitment(
        y = 1,
        rec_timeseries = r$full_recruitment,
        apportionment = r$rec_props,
        apportion_random = model_options$recruit_apportionment_random,
        apportionment_pars = model_options$recruit_apportionment_pars
    )

    expect_equal(r_y, array(20*c(0.5, 0.2, 0.1, 0.1, 0.1), dim=c(1, 5)))

})

test_that("Get annual recuits, fixed apportionment stochastic", {

    example_rec_generation_function <- function(input1, input2){
        return(c(0.5, 0.2, 0.1, 0.1, 0.1))
    }

    nyears <- 10
    nregions <- 5
    recruitment <- array(20, dim=c(nyears+1, 1))
    model_options <- list(
        recruit_apportionment = example_rec_generation_function,
        recruit_apportionment_pars = list(
            input1 = 12,
            input2 = 24
        ),
        recruit_apportionment_random = TRUE
    )

    r <- calculate_recruitment(recruitment, model_options$recruit_apportionment, nyears, nregions)

    set.seed(1120)
    r_y <- get_annual_recruitment(
        y = 1,
        rec_timeseries = r$full_recruitment,
        apportionment = r$rec_props,
        apportion_random = model_options$recruit_apportionment_random,
        apportionment_pars = model_options$recruit_apportionment_pars
    )

    expect_equal(r_y/sum(r_y), array(c(0.533, 0.133, 0.066, 0.133, 0.133), dim=c(1, 5)), tolerance=1e-2)
    expect_equal(r_y, array(c(10.67, 2.67, 1.33, 2.67, 2.67), dim=c(1, 5)), tolerance = 1e-1)

})
