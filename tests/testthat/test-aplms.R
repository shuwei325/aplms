library(testthat)

test_that("aplms halts the program when missing formula", {
  expect_error(
    aplms::aplms(npc=c("tdate","epi.week"), basis=c("cr","cc"),Knot=c(60,12),
                 data=hospitalization,family=Powerexp(k=0.3),p=3,
                 control = list(tol = 0.001,
                                algorithm1 = c("P-GAM"),
                                algorithm2 = c("BFGS"),
                                Maxiter1 = 20,
                                Maxiter2 = 25),
                 lam=c(100,10)),
    "The formula argument is missing."
  )
})

test_that("aplms halts the program when missing non parametric components", {
  expect_error(
    aplms::aplms(formula = y ~ 
                    MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max,
                  basis=c("cr","cc"),Knot=c(60,12),
                  data=hospitalization,family=Powerexp(k=0.3),p=3,
                  control = list(tol = 0.001,
                                 algorithm1 = c("P-GAM"),
                                 algorithm2 = c("BFGS"),
                                 Maxiter1 = 20,
                                 Maxiter2 = 25),
                  lam=c(100,10)),
    "The model needs at least one non-parametric component."
  )
})

test_that("aplms halts the program when missing data", {
  expect_error(
    aplms::aplms(formula = y ~ 
                    MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max,
                  npc=c("tdate","epi.week"), basis=c("cr","cc"),Knot=c(60,12),
                  family=Powerexp(k=0.3),p=3,
                  control = list(tol = 0.001,
                                 algorithm1 = c("P-GAM"),
                                 algorithm2 = c("BFGS"),
                                 Maxiter1 = 20,
                                 Maxiter2 = 25),
                  lam=c(100,10)),
    "The data argument is missing."
  )
})

test_that("the non-parametric variables must be in data", {
  expect_error(
    aplms::aplms(formula = y ~ 
                    MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max,
                  npc=c("mu","sigma"), basis=c("cr","cc"),Knot=c(60,12),
                  data=hospitalization,family=Powerexp(k=0.3),p=3,
                  control = list(tol = 0.001,
                                 algorithm1 = c("P-GAM"),
                                 algorithm2 = c("BFGS"),
                                 Maxiter1 = 20,
                                 Maxiter2 = 25),
                  lam=c(100,10)),
    "The non-parametric variables must be in data."
  )
})

test_that("Family must be from a known distribution", {
  bad_family <- list(a = 1)
  expect_error(
    aplms::aplms(formula = y ~ 
                    MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max,
                  npc=c("tdate","epi.week"), basis=c("cr","cc"),Knot=c(60,12),
                  data=hospitalization,family=bad_family,p=3,
                  control = list(tol = 0.001,
                                 algorithm1 = c("P-GAM"),
                                 algorithm2 = c("BFGS"),
                                 Maxiter1 = 20,
                                 Maxiter2 = 25),
                  lam=c(100,10)),
    "'family' not recognized"
  )
})

test_that("The funtion returns an object of class aplms", {
  data(hospitalization)
  mod2 <- aplms::aplms(formula = y ~ 
                    MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max,
                  npc=c("tdate","epi.week"), basis=c("cr","cc"),Knot=c(60,12),
                  data=hospitalization,family=Powerexp(k=0.3),p=3,
                  control = list(tol = 0.001,
                                 algorithm1 = c("P-GAM"),
                                 algorithm2 = c("BFGS"),
                                 Maxiter1 = 20,
                                 Maxiter2 = 25),
                  lam=c(100,10))
  expect_s3_class(mod2, "aplms")
})