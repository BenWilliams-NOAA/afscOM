test_that("no catch required", {

  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  tac <- 0.0
  F_f <- find_F(
      f_guess = 0.05, 
      naa     = naa,
      waa     = dem_params$waa,
      mort    = dem_params$mort,
      selex   = dem_params$sel[,,f],
      ret     = dem_params$ret[,,f],
      dmr     = dem_params$dmr[,,f],
      prov_catch = tac
  )

  expect_equal(F_f, 0.0, tolerance=1e-6)
})


test_that("Find F for TAC", {

  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  assessment <- dget("data/test.rdat")
  naa <- matrix(NA, nrow=model_params$nages, ncol=2)
  naa[,1] <- assessment$natage.female["2023",]*1e6
  naa[,2] <- assessment$natage.male["2023",]*1e6

  tac <- 4900
  suppressWarnings(
    F_f <- find_F(
      f_guess = 0.05, 
      naa     = naa,
      waa     = dem_params$waa,
      mort    = dem_params$mort,
      selex   = dem_params$sel[,,f],
      ret     = dem_params$ret[,,f],
      dmr     = dem_params$dmr[,,f],
      prov_catch = tac
  )
  )

  expect_equal(as.numeric(F_f), 9.566e-06, tolerance=1e-5)
})


test_that("Find F for TAC via bisection", {

  dem_params <- readRDS("data/sablefish_dem_matrices.RDS")
  model_params <- get_model_dimensions(dem_params$sel)
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  assessment <- dget("data/test.rdat")
  naa <- matrix(NA, nrow=model_params$nages, ncol=2)
  naa[,1] <- assessment$natage.female["2023",]*1e6
  naa[,2] <- assessment$natage.male["2023",]*1e6

  tac <- 4900
  suppressWarnings({
    F_f_bisections <- findF_bisection(
      f_guess = 0.05, 
      naa     = naa,
      waa     = dem_params$waa,
      mort    = dem_params$mort,
      selex   = dem_params$sel[,,f],
      ret     = dem_params$ret[,,f],
      dmr     = dem_params$dmr[,,f],
      prov_catch = tac
    ) 

    F_f_mle <- find_F(
      f_guess = 0.05, 
      naa     = naa,
      waa     = dem_params$waa,
      mort    = dem_params$mort,
      selex   = dem_params$sel[,,f],
      ret     = dem_params$ret[,,f],
      dmr     = dem_params$dmr[,,f],
      prov_catch = tac
    )

  })

  expect_equal(as.numeric(F_f_bisections), as.numeric(F_f_mle), tolerance=1e-5)
})


