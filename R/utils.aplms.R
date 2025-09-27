#' Calculate the first aproximation of the gamma parameters.
#'
#' @param k number of non parametric components.
#' @param y observed values of the response variable.
#' @param N_i function matrices.
#' @return initial aproximation of the gamma parameters.
calculatef_init <- function(k, y, N_i) {
  f_init <- vector("list", k + 1)
  f_init[[1]] <- rbind(mean(y), cbind(rep(0, dim(N_i[[1]])[2] - 1)))
  for (i in 1:k) {
    f_init[[i + 1]] <- cbind(rep(0, length = dim(N_i[[i + 1]])[2]))
  }
  return(f_init)
}

#' Estimate the variance the asymptotic covariance-variance matrix for the estimated
#' gamma parameters. 
#'
#' @param family probability distribution function.
#' @param phi vector of the estimated phi parameters.
#' @param const2 constant
#' @param k number of nn parametric components.
#' @param AN vector of the function matrices multiplied by the matrix A
#' @return estimation of the variance of f.
estimateVarF <- function(family, phi, const2, k, AN) {
  dg_t <- family$g2(args,
    df = family$df, r = family$r,
    s = family$s, alpha = family$alpha, mp = family$mp,
    epsi = family$epsi, sigmap = family$sigmap,
    k = family$k, nu = family$nu
  )

  const <- 4 * (dg_t / phi)
  rep_AN <- rep(AN, k + 1)
  seq_AN <- rep(AN, each = k + 1)

  fis_block <- mapply(FUN = function(X, Y) {
    t(X) %*% Y * const
  }, X = rep_AN, Y = seq_AN)

  diag_index <- diag(matrix(1:((k + 1)^2), ncol = k + 1))

  fis_block[diag_index] <- Map("+", fis_block[diag_index], const2)
  names(fis_block) <- letters[seq(from = 1, to = (k + 1)^2)]

  fis_FF <- list()
  for (l in 0:k) {
    fis_FF[[l + 1]] <- do.call(rbind, fis_block[(l * (k + 1) + 1):((l + 1) * (k + 1))])
  }

  fis_FF <- do.call(cbind, fis_FF)
  VAR_F <- solve(fis_FF)
  return(VAR_F)
}

#' Calculate the second derivative with respect to the gamma estimations
#'
#' @param phi vector of the estimated phi parameters.
#' @param Dd diagonal of the d1, ..., dn vectors
#' @param Dv diagonal of the v1, ..., vn vectors
#' @param AN vector of the function matrices multiplied by the matrix A
#' @param k number of non parametric components.
#' @param const2 constant
#' @return second derivative of the log function with respect to gamma
calculateLL_FF <- function(phi, Dd, Dv, AN, k, const2) {
  rep_AN <- rep(AN, k + 1)
  seq_AN <- rep(AN, each = k + 1)

  LL_block <- mapply(FUN = function(X, Y) {
    (1 / phi) * t(X) %*% (4 * Dd - Dv) %*% Y
  }, X = rep_AN, Y = seq_AN)

  diag_index <- diag(matrix(1:((k + 1)^2), ncol = k + 1))

  LL_block[diag_index] <- Map("-", LL_block[diag_index], const2)
  names(LL_block) <- letters[seq(from = 1, to = (k + 1)^2)]

  LL_FF <- list()
  for (l in 0:k) {
    LL_FF[[l + 1]] <- do.call(rbind, LL_block[(l * (k + 1) + 1):((l + 1) * (k + 1))])
  }

  LL_FF <- do.call(cbind, LL_FF)
  return(LL_FF)
}

#' Calculate the second derivative with respect to the gamma estimations and phi estimations.
#'
#' @param phi vector of the estimated phi parameters.
#' @param AN vector of the function matrices multiplied by the matrix A
#' @param Dd diagonal of the d1, ..., dn vectors
#' @param Dv diagonal of the v1, ..., vn vectors
#' @param A A matrix
#' @param error_hat residual of the model
#' @return second derivative of the log function with respect to tha gamma estimations and phi estimations
calculateLL_FF_Phi <- function(phi, AN, Dd, Dv, A, error_hat) {
  LL_FF_phi_block <- lapply(AN, FUN = function(x) {
    (1 / phi^2) * (t(x) %*% (2 * Dd - Dv) %*% (A %*% error_hat))
  })
  LL_FF_phi <- do.call(rbind, LL_FF_phi_block)
  return(LL_FF_phi)
}

#' Calculate the second derivative with respect to the phi estimations.
#'
#' @param nn number of observations.
#' @param phi vector of the estimated phi parameters.
#' @param delta_i delta_i
#' @param Dc diagonal of the  c1, ..., cn vectors
#' @param Dv diagonal of the  v1, ..., vn vectors
#' @return second derivative of the log function with respect to the phi estimations
calculateLL_Phi <- function(nn, phi, delta_i, Dc, Dv) {
  ONE <- cbind(rep(1, nn))
  LL_phi <- (1 / phi^2) * (nn / 2 + t(delta_i) %*% Dc %*% delta_i - t(delta_i) %*% Dv %*% ONE)
  return(LL_phi)
}

#' Calculate the second derivative with respect to the gamma estimations and rho estimations.
#' 
#' @param B B matrix
#' @param N_i function matrices.
#' @param p autoregressive order of the error
#' @param k number of non parametric components
#' @param phi vector of the estimated phi parameters.
#' @param A A matrix
#' @param Dd diagonal of the d1, ..., dn vectors
#' @param Dv diagonal of the  v1, ..., vn vectors
#' @param error_hat residual of the model
#' @return second derivative of the log function with respect to the gamma estimations and rho estimations
calculateLL_FF_Rho <- function(B, N_i, p, k, phi, A, Dv, Dd, error_hat) {
  rep_N_i <- rep(N_i, each = p)
  seq_B <- rep(B, k + 1)
  LL_FF_rho_block <- mapply(
    FUN = function(X, Y) {
      (1 / phi) * ((t(A %*% X) %*% (Dv - 4 * Dd) %*% (Y %*% error_hat)) + t(Y %*% X) %*% Dv %*% A %*% error_hat)
    },
    X = rep_N_i, Y = seq_B
  )
  LL_FF_rho <- list()
  for (l in 0:k) {
    LL_FF_rho[[l + 1]] <- do.call(cbind, LL_FF_rho_block[(l * (p) + 1):((l + 1) * p)])
  }
  LL_FF_rho <- do.call(rbind, LL_FF_rho)
  return(LL_FF_rho)
}

#' Calculate the second derivative with respect to the rho estimations.
#' 
#' @param B B matrix
#' @param p autoregressive order of the error
#' @param phi vector of the estimated phi parameters.
#' @param error_hat residual of the model
#' @param Dd diagonal of the d1, ..., dn vectors
#' @param Dv diagonal of the v1, ..., vn vectors
#' @return second derivative of the log function with respect to the rho estimations
calculateLL_Rho <- function(B, p, phi, error_hat, Dd, Dv) {
  rep_B <- rep(B, p)
  seq_B <- rep(B, each = p)
  LL_rho_block <- mapply(FUN = function(X, Y) {
    (1 / phi) * (t(X %*% error_hat) %*% (4 * Dd - Dv) %*% (Y %*% error_hat))
  }, X = rep_B, Y = seq_B)
  LL_rho <- matrix(LL_rho_block, nrow = p)
  return(LL_rho)
}

#' Calculate the second derivative with respect to the phi estimations and rho estimations.
#' 
#' @param B B matrix
#' @param phi vector of the estimated phi parameters.
#' @param error_hat residual of the model
#' @param Dd diagonal of the d1, ..., dn vectors
#' @param Dv diagonal of the v1, ..., vn vectors
#' @param A A matrix
#' @return second derivative of the log function with respect to the phi and rho estimations
calculateLL_Phi_Rho <- function(B, phi, error_hat, Dv, Dd, A) {
  LL_phi_rho_block <- lapply(B, FUN = function(x) {
    (1 / phi^2) * (t(x %*% error_hat) %*% (Dv - 2 * Dd) %*% (A %*% error_hat))
  })
  LL_phi_rho <- do.call(cbind, LL_phi_rho_block)
  return(LL_phi_rho)
}