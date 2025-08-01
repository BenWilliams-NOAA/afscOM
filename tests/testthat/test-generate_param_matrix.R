good_dem_matrices <- afscOM::good_dem_matrices
sablefish_assessment_data <- afscOM::sablefish_assessment_data

test_that("fill by single value", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )

  v <- 0.1
  m <- generate_param_matrix(v, dimension_names = dimension.names)
  expect_equal(m, good_dem_matrices$single)

})

test_that("fill by age", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )

  v <- exp(-5.1560+0.7331*2:31)/(1+exp(-5.1560+0.7331*2:31))
  m <- generate_param_matrix(v, dimension_names = dimension.names, by="age")
  expect_equal(m, good_dem_matrices$age)

})

test_that("fill by age and sex", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )

  log.waa.f <- log(5.87)+3.02*log(1-exp(-0.17*(2:31+2.98)))
  log.waa.m <- log(3.22)+3.02*log(1-exp(-0.27*(2:31+2.41)))
  waa.f <- exp(log.waa.f)
  waa.m <- exp(log.waa.m)

  v <- matrix(c(waa.f, waa.m), nrow=30)
  colnames(v) <- c("F", "M")
  rownames(v) <- 2:31

  m <- generate_param_matrix(v, dimension_names = dimension.names, by=c("age", "sex"))
  expect_equal(m, good_dem_matrices$age.sex)

})

test_that("fill by age, sex, and fleet", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )


  all_selex <- sablefish_assessment_data$agesel
  selex.f.ll_fish <- all_selex[,"fish1sel.f"]
  selex.m.ll_fish <- all_selex[,"fish1sel.m"]
  selex.f.tw_fish <- all_selex[,"fish3sel.f"]
  selex.m.tw_fish <- all_selex[,"fish3sel.m"]

  v <- array(NA, dim=c(30, 2, 2), dimnames=list(2:31, c("F", "M"), c("fixed", "trawl")))
  v[,1,1] <- selex.f.ll_fish
  v[,2,1] <- selex.m.ll_fish
  v[,1,2] <- selex.f.tw_fish
  v[,2,2] <- selex.m.tw_fish

  m <- generate_param_matrix(v, dimension_names = dimension.names, by=c("age", "sex", "fleet"))

  expect_equal(m, good_dem_matrices$age.sex.fleet, tolerance=1e-2)

})
