test_that("baranov catch equation is accurate", {
  fy <- 0.1
  naa <- array(c(100, 50, 25, 10), dim=c(1, 4, 1, 1))
  waa <- array(c(1, 2, 3, 4), dim=c(1, 4, 1, 1))
  mort <- array(rep(0.1, 4), dim=c(1, 4, 1, 1))
  selex <- array(c(0, 0.5, 1.0, 1.0), dim=c(1, 4, 1, 1, 1))

  suppressWarnings({
    catch <- baranov(fy, naa, waa, mort, selex)
  })

  expect_equal(catch, 15.066, tolerance=1e-4)
})
