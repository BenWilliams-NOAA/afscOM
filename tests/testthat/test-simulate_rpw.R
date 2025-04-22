assessment <- afscOM::sablefish_assessment_data
test_that("domestic longline survey RPW calculation", {


  ll_rpw_years <- as.numeric(rownames(assessment$obssrv1))
  ll_rpw_years <- ll_rpw_years-1960+1

  naa <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2, 1), dimnames = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"))
  naa[,,1,1] <- assessment$natage.female
  naa[,,2,1] <- assessment$natage.male

  ll_selex <- matrix(NA, nrow=30, ncol=2, dimnames=list(c(2:31), c("F", "M")))
  ll_selex[,1] <- assessment$agesel[,"srv1sel.f"]
  ll_selex[,2] <- assessment$agesel[,"srv1sel.m"]

  sel <- generate_param_matrix(ll_selex, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_survey"), by = c("age", "sex"), include_fleet_dim = TRUE)
  sel[(57:64),,1,,] <- matrix(rep(assessment$agesel[,"srv10sel.f"], length(27:34)), ncol=30, byrow=TRUE)
  sel[(57:64),,2,,] <- matrix(rep(assessment$agesel[,"srv10sel.m"], length(27:34)), ncol=30, byrow=TRUE)

  s <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2), dimnames=list("time"=1:nrow(assessment$natage.female), "age"=c(2:31), "sex"=c("F", "M")))
  s[,,1] <- assessment$surv.f
  s[,,2] <- assessment$surv.m
  survival <- generate_param_matrix(s, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"), by = c("time", "age", "sex"))
  zaa <- -log(survival)

  weight.mat <- as.matrix(assessment$growthmat[,c("wt.f.block1", "wt.m.block1")])
  dimnames(weight.mat) <- list("age"=c(2:31), "sex"=c("F", "M"))
  waa <- generate_param_matrix(weight.mat, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"), by=c("age", "sex"))

  ll_surv_q <- assessment$parameters["q1"]

  ll_preds <- assessment$obssrv1[,"predsrv1"]

  pred <- sapply(ll_rpw_years, function(y){
    simulate_rpw(ll_surv_q, naa[y,,,, drop=FALSE], waa[y,,,, drop=FALSE], subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5), zaa[y,,,, drop=FALSE])
  })

  pred/ll_preds

  expect_equal(as.vector(pred/ll_preds), rep(1, length(ll_rpw_years)), tolerance=1e-5)
})


test_that("GOA trawl survey RPW calculation", {
  tw_rpw_years <- c(1990, 1993, 1996, 1999, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023)-1960+1
  #tw_rpn_se    <- c(19.8844720848078, 18.9549441841442, 15.2997678164112, 22.2251374782261, 15.4044213762913, 16.964107652509, 17.3487199007805, 18.5446866076966, 20.3520955433346, 21.5679204610345, 23.1188125742554, 22.2790341441897, 29.1811814601308, 25.9353492635271, 27.1904485643506, 28.45526612239, 30.4573674832159, 25.5383376233216, 21.1452604831446, 25.2104581638705, 25.9218070073857, 30.5488014305134, 24.2567893977538, 20.4447712815498, 21.6188561193043, 18.5409375284709, 22.2018073450464, 31.7269054335811, 26.7913255523096, 51.4094728260081, 50.0462646144844, 50.0858647739984, 63.0175984493821, 74.8605699050248)

  naa <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2, 1), dimnames = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"))
  naa[,,1,1] <- assessment$natage.female
  naa[,,2,1] <- assessment$natage.male

  tw_selex <- matrix(NA, nrow=30, ncol=2, dimnames=list(c(2:31), c("F", "M")))
  tw_selex[,1] <- assessment$agesel[,"srv7sel.f"]
  tw_selex[,2] <- assessment$agesel[,"srv7sel.m"]

  sel <- generate_param_matrix(tw_selex, dimension_names = list("time"==1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="tw_survey"), by = c("age", "sex"), include_fleet_dim = TRUE)

  s <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2), dimnames=list("time"=1:nrow(assessment$natage.female), "age"=c(2:31), "sex"=c("F", "M")))
  s[,,1] <- assessment$surv.f
  s[,,2] <- assessment$surv.m
  survival <- generate_param_matrix(s, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"), by = c("time", "age", "sex"))
  zaa <- -log(survival)

  weight.mat <- as.matrix(assessment$growthmat[,c("wt.f.block1", "wt.m.block1")])
  dimnames(weight.mat) <- list("age"=c(2:31), "sex"=c("F", "M"))
  waa <- generate_param_matrix(weight.mat, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"), by=c("age", "sex"))

  tw_surv_q <- assessment$parameters["q7"]

  comp <- as.numeric(assessment$obssrv7[, "predsrv7"])

  pred <- sapply(tw_rpw_years, function(y){
    simulate_rpw(tw_surv_q, naa[y,,,, drop=FALSE], waa[y,,,, drop=FALSE], subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5), zaa[y,,,, drop=FALSE])
  })

  expect_equal(as.vector(pred/comp), rep(1, length(tw_rpw_years)), tolerance=1e-5)
})
