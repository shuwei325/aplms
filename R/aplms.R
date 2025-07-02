#' Fitting Additive partial linear models with symmetric errors
#'
#' \code{aplms} is used to fit additive partial linear models with symmetric errors.
#' In this setup, the natural cubic splines or cubic P-splines.
#'
#' @title Fitting Additive partial linear models with symmetric errors.
#' @param formula A symbolic description of the parametric component of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param npc A vector of names of non parametric component.
#' @param basis A vector of names of the basis to be used for each non parametric covariate.
#' @param Knot A vector of the number of knots in each non-linear component of the model.
#' @param data A data frame containing the variables in the model
#' @param family Symmetric error distribution
#' @param p autoregressive order of the error
#' @param control optimization rutine.
#' @param init A list of initial values for the symmetric error scale, phi, and autoregressive coefficients, rhos.
#' @param lam smoothing parameter vector.
#' @examples
#' a <- 1
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
  if (missingArg(formula)) {
    stop("The formula argument is missing.")
  }
  if (missingArg(npc)) {
    stop("The model needs at least one non-parametric component.")
  }
  if (missingArg(data)) {
    stop("The data argument is missing.")
  }
  if (!all(npc %in% names(data))) {
    stop("The non-parametric variables must be in data.")
  }
  if (control["algorithm2"] == "BFGS") {
    gradient <- NULL
  }
  # if (control["algorithm"]=="Fisher.score"){
  #   gradient=grr
  # }

  k <- length(npc)
  if (missingArg(basis)) {
    if (k == 1) basis <- c("ps")
    if (k >= 2) basis <- c("ps", rep("cp", k - 1)) # if (k==2) basis <- c("cr","cc")
  }

  if (missingArg(lam)) {
    lam <- rep(10, k)
  }

  data1 <- model.frame(formula, data = data)

  y <- model.response(data1)

  N0 <- model.matrix(formula, data = data1)
  q <- ncol(N0)
  nn <- nrow(N0)

  if (missingArg(Knot)) {
    Knot <- unlist(lapply(data[npc], max)) * (1 / 4)
    Knot <- sapply(Knot, function(x) min(x, 35))
  }

  K0 <- matrix(0, nrow = ncol(N0), ncol = ncol(N0))

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

  rdf <- nrow(data) - ncol(N0) - p - 1

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Symmetric error setup
  xi_t <- family$g4(1,
    df = family$df,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  )

  # initial values
  f_init <- vector("list", k + 1)
  f_init[[1]] <- rbind(mean(y), cbind(rep(0, dim(N_i[[1]])[2] - 1)))
  # assign nonparametric functions init.
  for (i in 1:k) {
    f_init[[i + 1]] <- cbind(rep(0, length = dim(N_i[[i + 1]])[2]))
  }

  if (missingArg(init)) {
    phi_ini <- sd(y) / xi_t
    rho_ini <- rep(0, p)
  } else {
    phi_ini <- init[[1]]
    rho_ini <- init[[2]]
  }

  f_aux <- f_init

  #####
  conv_betaf <- array()
  conv_betaf[1] <- 1
  conv_geral <- array()
  conv_geral[1] <- 1
  i <- j <- 1
  while (conv_geral[j] > control$tol) {
    print(paste("While", j))
    i <- 1
    conv_betaf <- array()
    conv_betaf[1] <- 1
    j <- j + 1
    A <- matrix_A(rho_ini, nn)


    if (control$algorithm1 == "backfitting") {
      while (conv_betaf[i] > control$tol) {
        a <- res(y, f_init, phi_ini, rho_ini, N_i)
        posicao <- as.vector(family$g1(a,
          df = family$df,
          alpha = family$alpha, mp = family$mp, epsi = family$epsi,
          sigmap = family$sigmap, k = family$k, nu = family$nu
        ))
        Dv <- (diag(-2 * posicao))
        ############################
        print(i)
        i <- i + 1
        if (i > control$Maxiter1) {
          break
        } ############################ control

        S_i <- list()
        S_i[[1]] <- tcrossprod(
          solve(t(A %*% N_i[[1]]) %*% Dv %*% (A %*% N_i[[1]])),
          (A %*% N_i[[1]])
        ) %*% Dv
        for (l in 1:k) {
          S_i[[l + 1]] <- tcrossprod(
            solve(t(A %*% N_i[[l + 1]]) %*% Dv %*% (A %*% N_i[[l + 1]]) + phi_ini * lam[l] * K_i[[l + 1]]),
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

        error <- mapply("-", f, f_init, SIMPLIFY = FALSE)
        conv_betaf[i] <- max(unlist(sapply(error, abs)))

        f_init <- f
      }
    } else if (control$algorithm1 == "P-GAM") {
      while (conv_betaf[i] > control$tol) {
        a <- res(y, f_init, phi_ini, rho_ini, N_i)
        posicao <- as.vector(family$g1(a,
          df = family$df,
          alpha = family$alpha, mp = family$mp, epsi = family$epsi,
          sigmap = family$sigmap, k = family$k, nu = family$nu
        ))
        Dv <- (diag(-2 * posicao))

        dg <- family$g2(a,
          df = family$df,
          alpha = family$alpha, mp = family$mp, epsi = family$epsi,
          sigmap = family$sigmap, k = family$k, nu = family$nu
        )
        ############################
        print(i)
        i <- i + 1
        if (i > control$Maxiter1) {
          break
        }

        AN_i <- lapply(N_i, FUN = function(x) {
          A %*% x
        })
        N_A <- do.call(cbind, AN_i)
        lam_list <- as.list(lam)
        lam_list <- append(0, lam_list)
        M_gamma <- mapply("*", lam_list, K_i, SIMPLIFY = FALSE)

        f_vector <- solve(t(N_A) %*% Dv %*% N_A + phi_ini * Matrix::bdiag(M_gamma)) %*% t(N_A) %*% Dv %*% (A %*% y)

        f <- list()
        param_f <- sapply(f_init, nrow)
        f_dimension <- cumsum(param_f)
        f[[1]] <- cbind(f_vector[1:f_dimension[1]])
        for (l in 2:length(f_dimension)) {
          f[[l]] <- cbind(f_vector[(f_dimension[l - 1] + 1):f_dimension[l]])
        }

        error <- mapply("-", f, f_init, SIMPLIFY = FALSE)
        conv_betaf[i] <- max(unlist(sapply(error, abs)))

        f_init <- f
      }
    }

    par1 <- optim(
      par = c(phi_ini, rho_ini),
      fn = logLik3.test, f = f_init, y = y, N_i = N_i, family = family,
      method = "L-BFGS-B",
      lower = c(0.001, rep(-1, p)), upper = c(Inf, rep(1, p)),
      control = list(fnscale = -1), hessian = T
    )

    dif_phi_rho <- par1$par - c(phi_ini, rho_ini)
    phi_ini <- par1$par[1]
    if (p > 0) {
      rho_ini <- par1$par[2:(1 + p)]
    }

    f_error <- mapply("-", f_aux, f, SIMPLIFY = FALSE)
    f_aux <- f

    conv_geral[j] <- max(abs(dif_phi_rho), unlist(sapply(f_error, abs)))

    if (j == control$Maxiter2) {
      break
    }
  }

  rho <- rho_ini
  phi <- phi_ini
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

  # loglik evaluation

  Lp <- logLik_fim.test(y, f, rho, phi, N_i, family)

  AN <- lapply(N_i, FUN = function(x) A %*% x)
  N_bar_a <- rlist::list.cbind(AN) # N_A
  N_bar <- rlist::list.cbind(N_i)
  K_ast <- phi * as.matrix(Matrix::bdiag(mapply("*", c(0, lam), K_i, SIMPLIFY = FALSE)))

  #
  H_alpha <- N_bar_a %*% solve((t(N_bar_a) %*% Dv) %*% N_bar_a + K_ast) %*% (t(N_bar_a) %*% Dv) # (4.18)

  q1 <- t(N_bar_a) %*% Dv %*% N_bar_a
  dec <- eigen(q1)
  q12 <- dec$vectors %*% (diag(dec$values^(1 / 2))) %*% t(dec$vectors)
  q12m <- dec$vectors %*% (diag(dec$values^(-1 / 2))) %*% t(dec$vectors)
  auto <- q12m %*% (K_ast) %*% q12m

  df_alpha <- sum(1 / (1 + eigen(auto)$value)) + p + 1 # Efective degree of freedom

  effect <- solve(t(N_bar_a) %*% Dv %*% N_bar_a + K_ast) %*% t(N_bar_a) %*% Dv %*% N_bar_a # pag 57, L4

  # effective degree of freedom per function
  n_i <- sapply(N_i, ncol)
  npc_dimension <- cumsum(n_i)
  dfk <- sum(diag(effect)[1:npc_dimension[1]])
  for (l in 2:length(npc_dimension)) {
    dfk_i <- sum(diag(effect)[(npc_dimension[l - 1] + 1):npc_dimension[l]])
    dfk <- append(dfk, dfk_i)
  }

  ######
  II <- diag(nn)
  yhat <- H_alpha %*% (A %*% y) - (A - II) %*% y

  AIC <- -2 * Lp + 2 * (df_alpha)
  BIC <- -2 * Lp + log(nn) * (df_alpha)
  AICC <- AIC + 1 +
    2 * (sum(diag(H_alpha)) + 1) / (nn - 2 - sum(diag(H_alpha)))
  AICc_alpha <- log((sum((sqrt(Dv) %*% (cbind(y - yhat)))^2)) / nn) + 1 +
    2 * (sum(diag(H_alpha)) + 1) / (nn - 2 - sum(diag(H_alpha)))


  muhat1 <- Reduce(`+`, mapply("%*%", N_i, f, SIMPLIFY = FALSE))
  error_hat <- y - muhat1

  GCV <- ((1 / nn) * sum(((sqrt(Dv) %*% (A %*% (y - muhat1)))^2))) /
    ((1 - sum(diag(H_alpha)) / nn)^2)



  ######################
  dg_t <- family$g2(args,
    df = family$df, r = family$r,
    s = family$s, alpha = family$alpha, mp = family$mp,
    epsi = family$epsi, sigmap = family$sigmap,
    k = family$k, nu = family$nu
  )
  fg_t <- family$g3(args,
    df = family$df, r = family$r,
    s = family$s, alpha = family$alpha, mp = family$mp,
    epsi = family$epsi, sigmap = family$sigmap,
    k = family$k, nu = family$nu
  )
  ###################################################

  Dd1 <- 4 * (dg_t / phi) * diag(nn)

  # Variancias de phi y rho.

  fis_phi <- (nn / (4 * phi^2)) * (4 * fg_t - 1) # pag 52
  # fis_rho = (4*dg_t*xi_t/(1-rho[1]^2))*(nn-1:p)

  VAR_phi <- 1 / fis_phi
  # VAR_rho<-1/fis_rho

  if (p > 0) {
    VAR_rho <- diag(solve(-par1$hessian))[2:(1 + p)]
  }

  # Variancias de las funciones suaves.

  const <- 4 * (dg_t / phi)
  const2 <- Map("*", K_i, c(0, lam))

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

  WALD <- list()
  for (l in 2:length(npc_dimension)) {
    WALD[[l - 1]] <- t(f[[l]]) %*% solve(VAR_F[(npc_dimension[l - 1] + 1):npc_dimension[l], (npc_dimension[l - 1] + 1):npc_dimension[l]]) %*% f[[l]]
  }

  WALD_vec <- unlist(WALD)

  WALD_p_value <- mapply(
    FUN = function(X, Y) {
      pchisq(q = X, df = Y, lower.tail = FALSE)
    },
    X = WALD_vec, Y = dfk[-1]
  )
  #

  ### LL_observada

  c_i <- as.vector(family$g5(a,
    df = family$df,
    alpha = family$alpha, mp = family$mp, epsi = family$epsi,
    sigmap = family$sigmap, k = family$k, nu = family$nu
  ))

  Dc <- diag(c_i)
  Dd <- diag(as.vector(a) * c_i)

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


  delta_i <- a^2
  ONE <- cbind(rep(1, nn))

  LL_phi <- (1 / phi^2) * (nn / 2 + t(delta_i) %*% Dc %*% delta_i - t(delta_i) %*% Dv %*% ONE)

  LL_FF_phi_block <- lapply(AN, FUN = function(x) {
    (1 / phi^2) * (t(x) %*% (2 * Dd - Dv) %*% (A %*% error_hat))
  })
  LL_FF_phi <- do.call(rbind, LL_FF_phi_block)

  if (p > 0) {
    B <- BB(p = p, nn = nn)
    rep_B <- rep(B, p)
    seq_B <- rep(B, each = p)
    LL_rho_block <- mapply(FUN = function(X, Y) {
      (1 / phi) * (t(X %*% error_hat) %*% (4 * Dd - Dv) %*% (Y %*% error_hat))
    }, X = rep_B, Y = seq_B)

    LL_rho <- matrix(LL_rho_block, nrow = p)

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


    LL_phi_rho_block <- lapply(B, FUN = function(x) {
      (1 / phi^2) * (t(x %*% error_hat) %*% (Dv - 2 * Dd) %*% (A %*% error_hat))
    })
    LL_phi_rho <- do.call(cbind, LL_phi_rho_block)

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

  # coeficientes

  f0 <- f[[1]]
  est_coef <- as.vector(f0)
  ee <- sqrt(diag(VAR_F)[1:length(f0)])
  t_test <- est_coef / ee

  p_value <- sapply(t_test, FUN = function(x) {
    2 * pt(abs(x), rdf, lower.tail = FALSE)
  })

  summary_table <- cbind(
    est_coef,
    ee,
    t_test,
    p_value
  )

  terms_formula <- stats::terms(formula)
  var_names <- attr(terms_formula, "term.labels")

  rownames(summary_table) <- c("intercept", var_names)
  colnames(summary_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  WALD_f <- cbind(WALD_vec, dfk[-1], WALD_p_value)
  rownames(WALD_f) <- npc
  colnames(WALD_f) <- c("Wald", "df", "Pr(>.)")

  if (p > 0) {
    WALD_rho <- rho / sqrt(VAR_rho)

    p_value_rho_normal <- sapply(WALD_rho, FUN = function(x) {
      2 * pnorm(abs(x), lower.tail = FALSE)
    })
    p_value_rho_t <- sapply(WALD_rho, FUN = function(x) {
      2 * pt(abs(x), rdf, lower.tail = FALSE)
    })
    summary_table_rho <- cbind(
      rho,
      sqrt(VAR_rho),
      WALD_rho,
      # p_value_rho_normal,
      p_value_rho_t
    )

    colnames(summary_table_rho) <- c("rho", "ee", "Wald", "p-value_t")
    rownames(summary_table_rho) <- c(paste0("rho", as.character(1:p)))
  } else {
    summary_table_rho <- NULL
  }

  WALD_phi <- phi / sqrt(VAR_phi)

  p_value_phi_normal <- 2 * pnorm(abs(WALD_phi), lower.tail = FALSE)
  p_value_phi_t <- 2 * pt(abs(WALD_phi), rdf, lower.tail = FALSE)
  summary_table_phi <- cbind(
    phi,
    sqrt(VAR_phi),
    WALD_phi,
    p_value_phi_t
  )

  rownames(summary_table_phi) <- "phi"


  summary_table_phirho <- rbind(summary_table_phi, summary_table_rho)
  colnames(summary_table_phirho) <- c("Estimate", "Std. Error", "Wald", "Pr(>|t|)")

  fit <- list(
    formula = formula, family = family, npc = npc, Knot = Knot,
    lam = lam, summary_table = summary_table, VAR_F = VAR_F,
    basis = basis,
    WALD_f = WALD_f,
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
    AICc_alpha = AICc_alpha,
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
#' @param model APLMS object.
#'
#' @param ... other arguments
#' @examples
#' a <- 1
#' @rdname aplms
#' @method print aplms
#' @export
print.aplms <- function(model, ...) {
  if (!inherits(model, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }
  cat("Call:\n")
  print(model$this.call)
}
