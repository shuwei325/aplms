#' Piece of the backfitting algorithm that creates the new estimation of the gamma parameters.
#' 
#' @param A A matrix
#' @param N_i function matrices
#' @param Dv diagonal of the v1, ..., vn vectors
#' @param k number of non parametric components
#' @param phi vector of the estimations of the phi parameters
#' @param lam smoothing parameter vector.
#' @param K_i function of matrices
#' @param f_init previous estimation of the gamma parameters
#' @param y vector of the observed values of the response variable
backfitting <- function(A, N_i, Dv, k, phi, lam, K_i, f_init, y) {
	S_i <- list()
	S_i[[1]] <- tcrossprod(
	  solve(t(A %*% N_i[[1]]) %*% Dv %*% (A %*% N_i[[1]])),
	  (A %*% N_i[[1]])
	) %*% Dv
	for (l in 1:k) {
	  S_i[[l + 1]] <- tcrossprod(
	    solve(t(A %*% N_i[[l + 1]]) %*% Dv %*% (A %*% N_i[[l + 1]]) + phi * lam[l] * K_i[[l + 1]]),
	    (A %*% N_i[[l + 1]])
	  ) %*% Dv
	}
	f <- f_init
	f0 <- S_i[[1]] %*% (A %*% (y - Reduce(`+`, mapply("%*%", N_i[-1], f[-1], SIMPLIFY = FALSE))))
	f[[1]] <- f0
	for (l in 1:k) {
	  f_i <- S_i[[l + 1]] %*% (A %*% (y - Reduce(`+`, mapply("%*%", N_i[-(l + 1)], f[-(l + 1)], SIMPLIFY = FALSE))))
	  f[[l + 1]] <- f_i
	}
	return(f)
}

#' Piece of the P-GAM algorithm that creates the new estimation of the gamma parameters.
#' 
#' @param A A matrix
#' @param N_i function matrices
#' @param Dv diagonal of the v1, ..., vn vectors
#' @param phi vector of the estimations of the phi parameters
#' @param lam smoothing parameter vector.
#' @param K_i function of matrices
#' @param f_init previous estimation of the gamma parameters
#' @param y vector of the observed values of the response variable
pgam <- function(A, N_i, Dv, phi, lam, K_i, f_init, y) {
	AN_i <- lapply(N_i, FUN = function(x) {
	  A %*% x
	})
	N_A <- do.call(cbind, AN_i)
	lam_list <- as.list(lam)
	lam_list <- append(0, lam_list)
	M_gamma <- mapply("*", lam_list, K_i, SIMPLIFY = FALSE)
	f_vector <- solve(t(N_A) %*% Dv %*% N_A + phi * Matrix::bdiag(M_gamma)) %*% t(N_A) %*% Dv %*% (A %*% y)
	f <- list()
	param_f <- sapply(f_init, nrow)
	f_dimension <- cumsum(param_f)
	f[[1]] <- cbind(f_vector[1:f_dimension[1]])
	for (l in 2:length(f_dimension)) {
	  f[[l]] <- cbind(f_vector[(f_dimension[l - 1] + 1):f_dimension[l]])
	}
	return(f)
}