test_that("domestic longline survey RPN calculation", {
  
  load(file.path(here::here(), "data/sablefish_assessment_data.rda"))
  assessment <- sablefish_assessment_data

  ll_rpn_years <- as.numeric(rownames(assessment$obssrv3))
  ll_rpn_years <- ll_rpn_years-1960+1

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

  ll_surv_q <- assessment$parameters["q1"]

  ll_preds <- assessment$obssrv3[,"predsrv3"]
  
  pred <- sapply(ll_rpn_years, function(y){
    simulate_rpn(ll_surv_q, naa[y,,,, drop=FALSE], subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5), zaa[y,,,, drop=FALSE])
  })

  pred/ll_preds

  expect_equal(as.vector(pred/ll_preds), rep(1, length(ll_rpn_years)), tolerance=1e-5)

})
