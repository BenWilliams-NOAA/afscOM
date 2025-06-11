dem_params <- afscOM::sablefish_dem_params

test_that("no catch required", {

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

  model_params <- get_model_dimensions(dem_params$sel)
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
  assessment <- sablefish_assessment_data
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

  model_params <- get_model_dimensions(dem_params$sel)
  y <- 1
  r <- 1
  f <- 1

  dem_params <- subset_dem_params(dem_params, y, d=1)
  dem_params <- subset_dem_params(dem_params, r, d=3)

  load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
  assessment <- sablefish_assessment_data
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

  expect_equal(as.numeric(F_f_bisections), as.numeric(F_f_mle), tolerance=1e-4)
})


test_that("Compare FAA to Sablefish Assessment", {

  assessment <- afscOM::sablefish_assessment_data

  TAC <- assessment$t.series[, "Catch_HAL"]

  naa <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2, 1), dimnames = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"))
  naa[,,1,1] <- assessment$natage.female
  naa[,,2,1] <- assessment$natage.male

  ll_selex <- matrix(NA, nrow=30, ncol=2, dimnames=list(c(2:31), c("F", "M")))
  ll_selex[,1] <- assessment$agesel[,"fish1sel.f"]
  ll_selex[,2] <- assessment$agesel[,"fish1sel.m"]
  sel <- generate_param_matrix(ll_selex, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_fishery"), by = c("age", "sex"), include_fleet_dim = TRUE)
  sel[(36:56),,1,,1] <- matrix(rep(assessment$agesel[, "fish4sel.f"], length(36:56)), ncol=30, byrow=TRUE)
  sel[(36:56),,2,,1] <- matrix(rep(assessment$agesel[, "fish4sel.m"], length(36:56)), ncol=30, byrow=TRUE)
  sel[(57:64),,1,,1] <- matrix(rep(assessment$agesel[, "fish5sel.f"], length(57:64)), ncol=30, byrow=TRUE)
  sel[(57:64),,2,,1] <- matrix(rep(assessment$agesel[, "fish5sel.m"], length(57:64)), ncol=30, byrow=TRUE)


  weight_mat <- matrix(NA, nrow=30, ncol=2, dimnames=list("age"=c(2:31), "sex"=c("F", "M")))
  weight_mat[,1] <- assessment$growthmat[, "wt.f.block1"]
  weight_mat[,2] <- assessment$growthmat[, "wt.m.block1"]
  waa <- generate_param_matrix(weight_mat, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_fishery"), by=c("age", "sex"))

  M <- 0.113179
  mort <- generate_param_matrix(M, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_fishery"))

  retention <- 1.0
  ret <- generate_param_matrix(retention, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_fishery"), include_fleet_dim = TRUE)

  discard <- 0.0
  dmr <- generate_param_matrix(discard, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_fishery"), include_fleet_dim = TRUE)

  fs  <- rep(NA, 64)
  faa <- array(NA, dim=c(64, 30, 2, 1))
  catch <- rep(NA, 64)

  fish1.fs <- c(0.00654294, 0.0343149, 0.0579996, 0.0238662, 0.00745146, 0.00201715, 0.00841816, 0.0086355, 0.0249528, 0.0352064, 0.0539089, 0.0575049, 0.0772958, 0.0684617, 0.081306, 0.0803652, 0.0956303, 0.076413, 0.0383202, 0.0435037, 0.0347983, 0.0442515, 0.0389907, 0.0359537, 0.0340402, 0.0402824, 0.0630109, 0.0799625, 0.0884755, 0.0895519, 0.0879386, 0.0829642, 0.0779795, 0.0889117, 0.0866495, 0.0714257, 0.0643549, 0.0592001, 0.0606085, 0.0603981, 0.0729324, 0.0670748, 0.0671666, 0.0747605, 0.0795493, 0.0720158, 0.0678445, 0.0710271, 0.0673023, 0.0634238, 0.0611044, 0.0691572, 0.0762706, 0.0764786, 0.0653089, 0.0637705, 0.0525657, 0.0555171, 0.0500482, 0.0438102, 0.0344363, 0.0401478, 0.04564, 0.037796)
  suppressWarnings({
    for(y in 1:64){
        F_f_bisections <- findF_bisection(
          f_guess = 0.05,
          naa     = naa[y,,,,drop=FALSE],
          waa     = waa[y,,,,drop=FALSE],
          mort    = mort[y,,,,drop=FALSE],
          selex   = subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5),
          ret     = subset_matrix(ret[y,,,,1, drop=FALSE], 1, d=5),
          dmr     = subset_matrix(dmr[y,,,,1, drop=FALSE], 1, d=5),
          prov_catch = TAC[y]
        )

        fs[y] <- F_f_bisections
        faa[y,,,] <- retained_F(F_f_bisections, subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5), subset_matrix(ret[y,,,,1, drop=FALSE], 1, d=5))
        catch[y] <- baranov(
          fy      = F_f_bisections,
          naa     = naa[y,,,,drop=FALSE],
          waa     = waa[y,,,,drop=FALSE],
          mort    = mort[y,,,,drop=FALSE],
          selex   = subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5),
          ret     = subset_matrix(ret[y,,,,1, drop=FALSE], 1, d=5),
          dmr     = subset_matrix(dmr[y,,,,1, drop=FALSE], 1, d=5)
        )
    }

  })
  dimnames(assessment$faa.fish1.f) <- list()

  expect_lt(max(faa[,,1,]-assessment$faa.fish1.f), 1e-4)
  expect_lt(max(fs - fish1.fs), 1e-4)
  expect_lt(max(catch - TAC), 1e-3)
})



