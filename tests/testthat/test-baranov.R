test_that("baranov catch equation is accurate", {
  fy <- 0.1
  naa <- c(100, 50, 25, 10)
  waa <- c(1, 2, 3, 4)
  mort <- rep(0.1, 4)
  selex <- c(0, 0.5, 1.0, 1.0)
  catch <- baranov(fy, naa, waa, mort, selex)

  expect_equal(catch, 15.066, tolerance=1e-4)
})
