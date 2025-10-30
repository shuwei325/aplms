# Utility Functions to support the likelihood function evaluation of the function `aplms`


# Compute the matrix A of autocorrelations parameters of the model.
#
# @param rho autoregressive parameters.
# @param nn dimension of the matrix A.
# @return the matrix A of autocorrelations.

matrix_A <- function(rho, nn) {
  if (length(rho) == 0) {
    A <- diag(nn)
  }
  if (length(rho) > 0) {
    p <- length(rho)
    A <- diag(nn)
    for (i in 2:(p + 1)) {
      for (j in 1:(nn + 1 - i)) {
        A[i + j - 1, j] <- -rho[i - 1]
      }
    }
  }
  return(A)
}

# Compute the matrix B of autocorrelations parameters of the model.
#
# @param P order of autoregressive parameters.
# @param nn dimension of the matrix B.
# @return the matrix B of autocorrelations.
BB <- function(p = 1, nn) {
  if (p == 0) {
    BB <- NULL
  }

  if (p > 0) {
    BB <- list()
    for (k in 1:p) {
      BB[[k]] <- matrix(0, nn, nn)
      for (i in (k + 1):nn) {
        BB[[k]][i, (i - k)] <- -1
      }
    }
  }
  return(BB)
}

# Residuals of the model to support the function aplms.
#
# @param resposta observed values of the response variable.
# @param eta estimated f.
# @param phi1 estimated scale parameter of the symmetric error.
# @param rho estimated autoregressive parameters of the model.
# @param N_i basis function matrices.
# @return residuals of the model.

res <- function(resposta, eta, phi1, rho, N_i) {
  nn <- length(resposta)
  A <- matrix_A(rho, nn)
  fitted1 <- mapply("%*%", N_i, eta, SIMPLIFY = FALSE)
  fitted2 <- Reduce(`+`, fitted1)
  val <- sqrt(1 / phi1) * A %*% (resposta - fitted2)
  val
}


# gradient second step
# grr<-function(par,f,y,N_i){ #
#   rho<-par[1:p]
#   phii1<-par[p+1]
#   A = matrix_A(rho,nn)
#   B = BB(p=p,nn)
#   one<-cbind(rep(1,nn))
#   epsi1<-epsi(y,  f_init,N_i)
#   a = res(y, f, par[1:p], par[p+1],N_i)
#   posicao = as.vector(family$g1(a, df = family$df,
#                                 alpha = family$alpha, mp = family$mp, epsi = family$epsi,
#                                 sigmap = family$sigmap, k = family$k))
#   Dv <- (diag(-2*posicao))
#   Dm <- diag(-2*posicao*c(a^2))
#   U <- numeric()
#   for(l in 1:p){
#     U[l] <- -1/phii1*t(B[[l]]%*%epsi1)%*%(Dv)%*%(A%*%epsi1)
#   }
#   return(c(U,t(one)%*%(Dm%*%one-one)/(2*phii1)))
# }


# profiled loglikelihood function evaluation to support the function aplms.
#
# @param par rhos and phi parameters to be optimized.
# @param f estimated f.
# @param y observed values of the response variable.
# @param N_i basis function matrices.
# @param family Symmetric error distribution.
# @return profiled loglikelihood function.

logLik3.test <- function(par, f, y, N_i, family) {
  p <- length(par) - 1
  nn <- length(y)
  phi <- par[1]
  rho <- par[-1]
  res <- res(y, f, phi, rho, N_i)
  logLik <- -(nn / 2) * log(par[1]) + sum(family$g0(res,
    df = family$df, s = family$s, r = family$r,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  ))
  logLik
}


# loglikelihood function evaluation to support the function aplms.
#
# @param y observed values of the response variable.
# @param f estimated f.
# @param rho estimated autoregressive parameters of the model.
# @param phi estimated scale parameter of the symmetric error.
# @param N_i basis function matrices.
# @param family Symmetric error distribution.
# @return estimated loglikelihood function.

logLik_fim.test <- function(y, f, rho, phi, N_i, family) {
  nn <- length(y)
  res <- res(y, f, phi, rho, N_i)
  logLik <- -(nn / 2) * log(phi) + sum(family$g0(res,
    df = family$df, s = family$s, r = family$r,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  ))
  logLik
}
