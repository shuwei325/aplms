#' @title Tests for the aplms utils function
#' @description 
#' Unit tests for the the aplms function. This tests include that the function returns
#' the correct answer when called with appropiate parameters and that the execution halts
#' when called with inappropiate parameters.
#' 

test_that("The return value of calculatef_init should be as expected", {
	data(hospitalization)
  expected_return <- list(as.matrix(c(83.50966, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000)),
                        as.matrix(rep(0, 59)), as.matrix(rep(0, 10)))
	
	formula <- y ~ MP10_avg + NO_avg + O3_avg + TEMP_min + ampl_max + RH_max
	npc <- c("tdate","epi.week")
	k <- length(npc)
	basis <- c("cr","cc") 
	Knot <- c(60,12)
  data <- hospitalization

	data1 <- model.frame(formula, data = data)
  y <- model.response(data1)
	N0 <- model.matrix(formula, data = data1)

  ZZ <- list()
  N_i <- list()
  K_i <- list()
  for (i in 1:k) {
    XX <- as.list(substitute(list(npc[i])))[-1]
    YY <- s(XX, bs = basis[i], m = c(2, 3), k = Knot[i])
    YY$term <- npc[i]
    ZZ[[i]] <- mgcv::smoothCon(YY, data = data, absorb.cons = T)
    N_i[[i]] <- ZZ[[i]][[1]]$X
  }

  N_i <- append(list(N0), N_i)

	f_init <- calculatef_init(k, y, N_i)
	
  expect_equal(expected_return, f_init, tolerance = 1e-6)
})

test_that("The return value of VAR_F should be as expected", {
	family <- readRDS("data_tests/parameters/parameter_family.rds")
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	const2 <- readRDS("data_tests/parameters/parameter_const2.rds")
	AN <- readRDS("data_tests/parameters/parameter_AN.rds")
	k <- 2

	VAR_F <- estimateVarF(family, phi, const2, k, AN)

	expected_return <- readRDS("data_tests/expected_values/expected_value_VAR_F.rds")

  expect_equal(expected_return, VAR_F, tolerance = 1e-6)
})

test_that("The return value of LL_FF should be as expected", {
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dd <- readRDS("data_tests/parameters/parameter_Dd.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	AN <- readRDS("data_tests/parameters/parameter_AN.rds")
	const2 <- readRDS("data_tests/parameters/parameter_const2.rds")
	k <- 2

	LL_FF <- calculateLL_FF(phi, Dd, Dv, AN, k, const2)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_FF.rds")

  expect_equal(expected_return, LL_FF, tolerance = 1e-6)
})

test_that("The return value of LL_FF_Phi should be as expected", {
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dd <- readRDS("data_tests/parameters/parameter_Dd.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	AN <- readRDS("data_tests/parameters/parameter_AN.rds")
	A <- readRDS("data_tests/parameters/parameter_A.rds")
	error_hat <- readRDS("data_tests/parameters/parameter_error_hat.rds")

	LL_FF_Phi <-calculateLL_FF_Phi(phi, AN, Dd, Dv, A, error_hat)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_FF_phi.rds")

  expect_equal(expected_return, LL_FF_Phi, tolerance = 1e-6)
})

test_that("The return value of LL_Phi should be as expected", {
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	Dc <- readRDS("data_tests/parameters/parameter_Dc.rds")
	delta_i <- readRDS("data_tests/parameters/parameter_delta_i.rds")
	nn <- readRDS("data_tests/parameters/parameter_nn.rds")

	LL_Phi <- calculateLL_Phi(nn, phi, delta_i, Dc, Dv)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_phi.rds")

  expect_equal(expected_return, LL_Phi, tolerance = 1e-6)
})

test_that("The return value of LL_FF_Rho should be as expected", {
	B <- readRDS("data_tests/parameters/parameter_B.rds")
	N_i <- readRDS("data_tests/parameters/parameter_N_i.rds")
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dd <- readRDS("data_tests/parameters/parameter_Dd.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	A <- readRDS("data_tests/parameters/parameter_A.rds")
	error_hat <- readRDS("data_tests/parameters/parameter_error_hat.rds")
	k <- 2
	p <- 3

	LL_FF_Rho <- calculateLL_FF_Rho(B, N_i, p, k, phi, A, Dv, Dd, error_hat)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_FF_rho.rds")

  expect_equal(expected_return, LL_FF_Rho, tolerance = 1e-6)
})

test_that("The return value of LL_Rho should be as expected", {
	B <- readRDS("data_tests/parameters/parameter_B.rds")
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dd <- readRDS("data_tests/parameters/parameter_Dd.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	error_hat <- readRDS("data_tests/parameters/parameter_error_hat.rds")
	p <- 3

	LL_Rho <- calculateLL_Rho(B, p, phi, error_hat, Dd, Dv)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_rho.rds")
  expect_equal(expected_return, LL_Rho, tolerance = 1e-6)
})

test_that("The return value of LL_Phi_Rho should be as expected", {
	B <- readRDS("data_tests/parameters/parameter_B.rds")
	phi <- readRDS("data_tests/parameters/parameter_phi.rds")
	Dd <- readRDS("data_tests/parameters/parameter_Dd.rds")
	Dv <- readRDS("data_tests/parameters/parameter_Dv.rds")
	A <- readRDS("data_tests/parameters/parameter_A.rds")
	error_hat <- readRDS("data_tests/parameters/parameter_error_hat.rds")

	LL_phi_rho <- calculateLL_Phi_Rho(B, phi, error_hat, Dv, Dd, A)

	expected_return <- readRDS("data_tests/expected_values/expected_value_LL_phi_rho.rds")

  expect_equal(expected_return, LL_phi_rho, tolerance = 1e-6)
})