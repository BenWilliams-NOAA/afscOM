test_that("fill by single value", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )

  v <- 0.1
  m <- generate_param_matrix(v, dimension.names = dimension.names)
  good <- readRDS("data/good_dem_matrices.RDS")

  expect_equal(m, good$single)

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
  m <- generate_param_matrix(v, dimension.names = dimension.names, by="age")
  good <- readRDS("data/good_dem_matrices.RDS")

  expect_equal(m, good$age)

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
  rownames(waa.matrix) <- 2:31

  m <- generate_param_matrix(waa.matrix, dimension.names = dimension.names, by=c("age", "sex")) 
  good <- readRDS("data/good_dem_matrices.RDS")

  expect_equal(m, good$age.sex)

})

test_that("fill by age, sex, and fleet", {

  dimension.names <- list(
    "time" = 1:10,
    "age" = 2:31,
    "sex" = c("F", "M"),
    "region" = c("GOA", "BSAI"),
    "fleet" = c("fixed", "trawl")
  )

  all_selex <- dget("data/test.rdat")$agesel
  selex.f.ll_fish <- all_selex[,"fish1sel.f"]
  selex.m.ll_fish <- all_selex[,"fish1sel.m"]
  selex.f.tw_fish <- all_selex[,"fish3sel.f"]
  selex.m.tw_fish <- all_selex[,"fish3sel.m"]

  v <- array(NA, dim=c(30, 2, 2), dimnames=list(2:31, c("F", "M"), c("fixed", "trawl")))
  v[,1,1] <- selex.f.ll_fish
  v[,2,1] <- selex.m.ll_fish
  v[,1,2] <- selex.f.tw_fish
  v[,2,2] <- selex.m.tw_fish

  m <- generate_param_matrix(v, dimension.names = dimension.names, by=c("age", "sex", "fleet"))
  good <- readRDS("data/good_dem_matrices.RDS")

  expect_equal(m, good$age.sex.fleet)

})
