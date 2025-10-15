#' @title Fitting Additive partial linear models with symmetric errors
#' @description
#' \code{aplms} fits additive partial linear models with autoregressive symmetric errors.
#' This method is suitable for data sets where the response variable is continuous and symmetric,
#' with either heavy or light tails, and measured over time.
#' The model includes a parametric component for a set of covariates, while another set of covariates
#' can be specified as semi-parametric functions, typically time-related.
#' In this setup, natural cubic splines or cubic P-splines are used to approximate the nonparametric components.
#' @param formula A symbolic description of the parametric component of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param npc A vector of names of non parametric component.
#' @param basis A vector of names of the basis to be used for each non parametric covariate.
#' @param Knot A vector of the number of knots in each non-linear component of the model.
#' @param data A data frame containing the variables in the model
#' @param family Symmetric error distribution. The implemented distribution are: \code{Normal()}, \code{LogisI()}, \code{LogisII()}, \code{Student(df)}, \code{Powerexp(k)}, \code{Gstudent(parm=c(s,r))}.
#' @param p autoregressive order of the error
#' @param control optimization rutine.
#' @param init A list of initial values for the symmetric error scale, phi, and autoregressive coefficients, rhos.
#' @param lam smoothing parameter vector.
#' @return Returns an object of class \dQuote{aplms}, a list with following components.
#' \item{formula}{the \code{formula} object used.}
#' \item{family}{the \code{family} object used.}
#' \item{npc}{the \code{npc} object used.}
#' \item{Knot}{the \code{Knot} object used.}
#' \item{lam}{the \code{lam} object used.}
#' \item{rdf}{Degrees of freedom: \code{n - q - p - 1}.}
#' \item{VAR_F}{Estimate the asymptotic covariance matrix for the gamma parameters.}
#' \item{basis}{The \code{basis} to be used for each non parametric covariate.}
#' \item{WALD_f}{The summary table of the Wald statistics.}
#' \item{summary_table_phirho}{The summary table of the rho and phi parameters.}
#' \item{N_i}{Basis functions.}
#' \item{f}{Estimated gamma parameters.}
#' \item{Dv}{Dv values for the symmetric error.}
#' \item{Dm}{Dm values for the symmetric error.}
#' \item{Dc}{Dc values for the symmetric error.}
#' \item{Dd}{Dd values for the symmetric error.}
#' \item{delta}{delta_i for the symmetric error.}
#' \item{LL_obs}{Observed information matrix of the fitted model.}
#' \item{loglike}{The estimated loglikelihood function of the fitted model.}
#' \item{total_df}{The total effective degree of freedom of the model.}
#' \item{parametric_df}{The degree of freedom of the parametric components.}
#' \item{npc_df}{The effective degree of freedom of the non parametric components.}
#' \item{AIC}{Akaike information criterion of the estimated model.}
#' \item{BIC}{Bayesian information criterion of the estimated model.}
#' \item{AICC}{Corrected Akaike information criterion of the estimated model.}
#' \item{GCV}{The generalized cross-validation (GCV).}
#' \item{yhat}{The fitted response values of the model.}
#' \item{muhat}{The fitted mean values of the model.}
#' \item{residuals_y}{The response residuals}
#' \item{residuals_mu}{Raw (Ordinary) residuals: \eqn{y_t - (\textbf{x}_i^\top\beta + f_1(t_{i1}) + \ldots + f_k(t_{ik}))}}
#' \item{data}{the \code{data} object used.}
#' \item{this.call}{the function call used.}
#' @examples
#' data(temperature)
#' datos = data.frame(temperature,time=1:length(temperature))
#' mod1<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=datos,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' summary(mod1)
#' print(mod1)
#' @references Chou-Chen, S.W., Oliveira, R.A., Raicher, I., Gilberto A. Paula (2024)
#'    \emph{Additive partial linear models with autoregressive symmetric errors and its application to the
#'    hospitalizations for respiratory diseases} Stat Papers 65, 5145–5166. \doi{10.1007/s00362-024-01590-w}
#' @references Oliveira, R.A., Paula, G.A. (2021)
#'    \emph{Additive partial linear models with autoregressive symmetric errors and its application to the
#'    hospitalizations for respiratory diseases} Comput Stat 36, 2435–2466. \doi{10.1007/s00180-021-01106-2}
#' @encoding UTF-8
#' @import methods mgcv stats
#' @export aplms
aplms <- function(formula, npc, basis, Knot, data, family = Normal(), p = 1,
                  control = list(
                    tol = 0.001,
                    algorithm1 = c("P-GAM"),
                    algorithm2 = c("BFGS"),
                    Maxiter1 = 20,
                    Maxiter2 = 25
                  ),
                  init,
                  lam) {
  this.call <- match.call()
  if (missingArg(formula)) stop("The formula argument is missing.")
  if (missingArg(npc)) stop("The model needs at least one non-parametric component.")
  if (missingArg(data)) stop("The data argument is missing.")
  if (!all(npc %in% names(data))) stop("The non-parametric variables must be in data.")
  if (is.null(family$family)) stop("'family' not recognized")
  if (control$algorithm1 != "P-GAM" & control$algorithm1 != "backfitting") stop("The algorithm should be P-GAM or backfitting.")
  if(control$Maxiter1 < 0 | control$Maxiter2 < 0) stop("Maxiter1 and Maxiter2 should be positive integers.")
  if(control$tol < 0) stop("The tolerance should be positive number.")

  k <- length(npc)

  if (missingArg(basis)) {
    if (k == 1) basis <- c("ps")
    if (k >= 2) basis <- c("ps", rep("cp", k - 1))
  } else {
    if (length(basis) != k) stop("The vector of names of the basis should be the same length as the non-parametric component.")
  }

  if (missingArg(lam)) {
    lam <- rep(10, k)
  } else {
    if (length(lam) != k) stop("The smoothing parameter vector should be the same length as the non-parametric component.")
  }

  if (missingArg(Knot)) {
    Knot <- unlist(lapply(data[npc], max)) * (1 / 4)
    Knot <- sapply(Knot, function(x) min(x, 35))
  } else {
    if (length(Knot) != k) stop("The vector of the knots should be the same length as the non-parametric component.")
  }

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  # Symmetric error setup
  xi_t <- family$g4(1,
    df = family$df,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  )
  data1 <- model.frame(formula, data = data)

  y <- model.response(data1)
  if (missingArg(init)) {
    phi <- sd(y) / xi_t
    rho <- rep(0, p)
  } else {
    if(length(init) != 2) stop("There should be an initial value for the symmetric error scale and for the autoregressive coefficients")
    phi <- init[[1]]
    rho <- init[[2]]
  }

  N0 <- model.matrix(formula, data = data1)
  q <- ncol(N0)
  nn <- nrow(N0)


  K0 <- matrix(0, nrow = q, ncol = q)

  ZZ <- list()
  N_i <- list()
  K_i <- list()
  for (i in 1:k) {
    XX <- as.list(substitute(list(npc[i])))[-1]
    YY <- s(XX, bs = basis[i], m = c(2, 3), k = Knot[i])
    YY$term <- npc[i]
    ZZ[[i]] <- mgcv::smoothCon(YY, data = data, absorb.cons = T)
    N_i[[i]] <- ZZ[[i]][[1]]$X
    K_i[[i]] <- (ZZ[[i]][[1]]$S)[[1]]
  }

  N_i <- append(list(N0), N_i)
  K_i <- append(list(K0), K_i)

  rdf <- nrow(data) - q - p - 1
  f_init <- calculatef_init(k, y, N_i)
  f_aux <- f_init

  #####
  conv_geral <- 1
  j <- 1
  while (conv_geral > control$tol && j < control$Maxiter2) {
    print(paste("Iteration", j ))
    i <- 1
    conv_betaf <- 1
    A <- matrix_A(rho, nn)

    while (conv_betaf > control$tol && i < control$Maxiter1) {
      cat(paste("Iteration", j,"-", i))
      print(i)
      a <- res(y, f_init, phi, rho, N_i)
      posicao <- as.vector(family$g1(a,
        df = family$df,
        alpha = family$alpha, mp = family$mp, epsi = family$epsi,
        sigmap = family$sigmap, k = family$k, nu = family$nu
      ))
      Dv <- (diag(-2 * posicao))

      if (control$algorithm1 == "backfitting") {
        f <- backfitting(A, N_i, Dv, k, phi, lam, K_i, f_init, y)
      } else if (control$algorithm1 == "P-GAM") {
        f <- pgam(A, N_i, Dv, phi, lam, K_i, f_init, y)
      }

      error <- mapply("-", f, f_init, SIMPLIFY = FALSE)
      conv_betaf <- max(unlist(sapply(error, abs)))
      f_init <- f
      i <- i + 1
    }

    par1 <- optim(
      par = c(phi, rho),
      fn = logLik3.test, f = f_init, y = y, N_i = N_i, family = family,
      method = "L-BFGS-B",
      lower = c(0.001, rep(-1, p)), upper = c(Inf, rep(1, p)),
      control = list(fnscale = -1), hessian = T
    )

    dif_phi_rho <- par1$par - c(phi, rho)
    phi <- par1$par[1]
    if (p > 0) {
      rho <- par1$par[2:(1 + p)]
    }

    f_error <- mapply("-", f_aux, f, SIMPLIFY = FALSE)
    f_aux <- f

    conv_geral <- max(abs(dif_phi_rho), unlist(sapply(f_error, abs)))
    j <- j + 1
  }

  A <- matrix_A(rho, nn)
  a <- res(y, f, phi, rho, N_i)
  posicao <- as.vector(family$g1(a,
    df = family$df,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  ))
  Dv <- (diag(-2 * posicao))
  Dm <- diag(-2 * posicao * as.vector(a)^2)
  Dd <- diag(posicao * as.vector(a))
  delta_i <- a^2

  # loglik evaluation
  Lp <- logLik_fim.test(y, f, rho, phi, N_i, family)
  AN <- lapply(N_i, FUN = function(x) A %*% x)
  N_bar_a <- rlist::list.cbind(AN)
  K_ast <- phi * as.matrix(Matrix::bdiag(mapply("*", c(0, lam), K_i, SIMPLIFY = FALSE)))

  # effective degree of freedom per function
  effect <- solve(t(N_bar_a) %*% Dv %*% N_bar_a + K_ast) %*% t(N_bar_a) %*% Dv %*% N_bar_a # pag 57, L4
  n_i <- sapply(N_i, ncol)
  npc_dimension <- cumsum(n_i)
  dfk <- sum(diag(effect)[1:npc_dimension[1]])
  for (l in 2:length(npc_dimension)) {
    dfk_i <- sum(diag(effect)[(npc_dimension[l - 1] + 1):npc_dimension[l]])
    dfk <- append(dfk, dfk_i)
  }

  ######
  q1 <- t(N_bar_a) %*% Dv %*% N_bar_a
  dec <- eigen(q1)
  q12 <- dec$vectors %*% (diag(dec$values^(1 / 2))) %*% t(dec$vectors)
  q12m <- dec$vectors %*% (diag(dec$values^(-1 / 2))) %*% t(dec$vectors)
  auto <- q12m %*% (K_ast) %*% q12m
  df_alpha <- sum(1 / (1 + eigen(auto)$value)) + p + 1 # Efective degree of freedom

  II <- diag(nn)
  H_alpha <- N_bar_a %*% solve((t(N_bar_a) %*% Dv) %*% N_bar_a + K_ast) %*% (t(N_bar_a) %*% Dv) # (4.18)
  yhat <- H_alpha %*% (A %*% y) - (A - II) %*% y

  AIC <- -2 * Lp + 2 * (df_alpha)
  BIC <- -2 * Lp + log(nn) * (df_alpha)
  AICC <- AIC + 1 +
    2 * (sum(diag(H_alpha)) + 1) / (nn - 2 - sum(diag(H_alpha)))

  muhat1 <- Reduce(`+`, mapply("%*%", N_i, f, SIMPLIFY = FALSE))
  error_hat <- y - muhat1

  GCV <- ((1 / nn) * sum(((sqrt(Dv) %*% (A %*% (y - muhat1)))^2))) /
    ((1 - sum(diag(H_alpha)) / nn)^2)

  ### LL_observada
  c_i <- as.vector(family$g5(a,
    df = family$df,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  ))

  Dc <- diag(c_i)
  Dd <- diag(as.vector(a) * c_i)
  const2 <- Map("*", K_i, c(0, lam))
  LL_FF <- calculateLL_FF(phi, Dd, Dv, AN, k, const2)
  LL_FF_phi <- calculateLL_FF_Phi(phi, AN, Dd, Dv, A, error_hat)
  LL_phi <- calculateLL_Phi(nn, phi, delta_i, Dc, Dv)

  if (p > 0) {
    B <- BB(p = p, nn = nn)
    LL_rho <- calculateLL_Rho(B, p, phi, error_hat, Dd, Dv)
    LL_FF_rho <- calculateLL_FF_Rho(B, N_i, p, k, phi, A, Dv, Dd, error_hat)
    LL_phi_rho <- calculateLL_Phi_Rho(B, phi, error_hat, Dv, Dd, A)
    LL_obs <- rbind(
      cbind(LL_FF, LL_FF_phi, LL_FF_rho),
      cbind(t(LL_FF_phi), LL_phi, LL_phi_rho),
      cbind(t(LL_FF_rho), t(LL_phi_rho), LL_rho)
    )
  } else {
    LL_obs <- rbind(
      cbind(LL_FF, LL_FF_phi),
      cbind(t(LL_FF_phi), LL_phi)
    )
  }

  # Summary Tables
  VAR_F <- estimateVarF(family, phi, const2, k, AN)
  WALD_f <- generateWaldF(npc_dimension, dfk, npc, f, VAR_F)
  summary_table_phirho <- generateSummartTablePhiRho(p, par1, rdf, rho, phi, nn, family)

  fit <- list(
    formula = formula, family = family, npc = npc, Knot = Knot,
    lam = lam, rdf = rdf, VAR_F = VAR_F, basis = basis, WALD_f = WALD_f,
    summary_table_phirho = summary_table_phirho,
    N_i = N_i, f = f,
    Dv = Dv,
    Dm = Dm,
    Dc = Dc,
    Dd = Dd,
    delta = delta_i,
    LL_obs = LL_obs,
    loglike = Lp,
    total_df = df_alpha,
    parametric_df = dfk[1],
    npc_df = dfk[-1],
    AIC = AIC,
    BIC = BIC,
    AICC = AICC,
    #AICc_alpha = AICc_alpha,
    GCV = GCV,
    yhat = yhat,
    muhat = muhat1,
    residuals_y = y - yhat,
    residuals_mu = y - muhat1,
    data = data,
    this.call = this.call
  )
  class(fit) <- "aplms"
  return(fit)
}

#' Print APLMS object.
#' @param x APLMS object.
#' @param ... other arguments
#' @rdname aplms
#' @export
print.aplms <- function(x, ...) {
  if (!inherits(x, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }
  cat("Call:\n")
  print(x$this.call)
}
