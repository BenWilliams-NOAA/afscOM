test_that("domestic longline survey RPN calculation", {
  
  ll_rpn_years <- c(1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)-1960+1
  ll_rpn_se    <- c(19.8844720848078, 18.9549441841442, 15.2997678164112, 22.2251374782261, 15.4044213762913, 16.964107652509, 17.3487199007805, 18.5446866076966, 20.3520955433346, 21.5679204610345, 23.1188125742554, 22.2790341441897, 29.1811814601308, 25.9353492635271, 27.1904485643506, 28.45526612239, 30.4573674832159, 25.5383376233216, 21.1452604831446, 25.2104581638705, 25.9218070073857, 30.5488014305134, 24.2567893977538, 20.4447712815498, 21.6188561193043, 18.5409375284709, 22.2018073450464, 31.7269054335811, 26.7913255523096, 51.4094728260081, 50.0462646144844, 50.0858647739984, 63.0175984493821, 74.8605699050248)

  assessment <- dget("data/test.rdat")
  naa <- array(NA, dim=c(length(ll_rpn_years), 30, 2, 1), dimnames = list("time"=ll_rpn_years, "age"=2:31, "sex"=c("F", "M"), "region"="alaska"))
  naa[,,1,1] <- assessment$natage.female[ll_rpn_years,]
  naa[,,2,1] <- assessment$natage.male[ll_rpn_years,]

  ll_selex <- matrix(NA, nrow=30, ncol=2, dimnames=list(c(2:31), c("F", "M")))
  ll_selex[,1] <- assessment$agesel[,"srv1sel.f"]
  #ll_selex[,1] <- c(0.00442, 0.06400, 0.51200, 0.94151, 0.99596, 0.99973, 0.99998, 0.99999, rep(1, 22))
  ll_selex[,2] <- assessment$agesel[,"srv1sel.m"]
  #ll_selex[,2] <- c(0.01879, 0.17728, 0.70796, 0.96463, 0.99675, 0.99971, 0.99997, 0.99999, rep(1, 22))

  sel <- generate_param_matrix(ll_selex, dimension_names = list("time"=ll_rpn_years, "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_survey"), by = c("age", "sex"), include_fleet_dim = TRUE)
  sel[(57:64)-30, , 1,, ] <- assessment$agesel[, "srv10sel.f"]
  #sel[(57:64)-30,,1,,] <- matrix(rep(c(0.08147, 0.57735, 0.9544, 0.9969, 0.99979, 0.99998, 0.99999, rep(1, 23)), length(57:64)), ncol=30, byrow=TRUE)
  sel[(57:64)-30, , 2,, ] <- assessment$agesel[, "srv10sel.m"]
  #sel[(57:64)-30,,2,,] <- matrix(rep(c(0.14554, 0.65711, 0.9556, 0.9958, 0.99963, 0.99996, 0.99999, rep(1, 23)), length(57:64)), ncol=30, byrow=TRUE)

  ll_surv_q <- 6.41/1.09

  comp <- as.numeric(assessment$obssrv3[, "predsrv3"])
  
  pred <- sapply(ll_rpn_years-30, function(y){
    simulate_rpn(ll_surv_q, naa[y,,,, drop=FALSE], subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5))
  })

  expect_equal(pred/comp, rep(1, length(ll_rpn_years)), tolerance=0.1)

})
