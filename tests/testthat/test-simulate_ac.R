test_that("domestic LL survey age composition calculation", {
  
  assessment <- dget("data/test.rdat")

  ll_rpn_years <- as.numeric(rownames(assessment$eac.srv1))
  ll_rpn_years <- ll_rpn_years-1960+1

  naa <- array(NA, dim=c(nrow(assessment$natage.female), 30, 2, 1), dimnames = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska"))
  naa[,,1,1] <- assessment$natage.female
  naa[,,2,1] <- assessment$natage.male

  ll_selex <- matrix(NA, nrow=30, ncol=2, dimnames=list(c(2:31), c("F", "M")))
  ll_selex[,1] <- assessment$agesel[,"srv1sel.f"]
  ll_selex[,2] <- assessment$agesel[,"srv1sel.m"]

  sel <- generate_param_matrix(ll_selex, dimension_names = list("time"=1:nrow(assessment$natage.female), "age"=2:31, "sex"=c("F", "M"), "region"="alaska", "fleet"="ll_survey"), by = c("age", "sex"), include_fleet_dim = TRUE)
  sel[(57:64),,1,,] <- matrix(rep(assessment$agesel[,"srv10sel.f"], length(57:64)), ncol=30, byrow=TRUE)
  sel[(57:64),,2,,] <- matrix(rep(assessment$agesel[,"srv10sel.m"], length(57:64)), ncol=30, byrow=TRUE)

  ageage <- read_csv("~/Desktop/ageage.csv", show_col_types=FALSE) %>% column_to_rownames("...1") %>% as.matrix

  ll_preds <- assessment$eac.srv1
  dimnames(ll_preds) <- list()

  pred <- sapply(ll_rpn_years, function(y){
    simulate_ac(naa[y,,,, drop=FALSE], subset_matrix(sel[y,,,,1, drop=FALSE], 1, d=5), ageage)
  })

  expect_equal(t(pred), ll_preds, tolerance=1e-5)
})
