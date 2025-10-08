
#' Print method for "aplms" class
#'
#' @param object an object with the result of fitting additive partial linear models with symmetric errors.
#' @param ... Other arguments passed to or from other methods.
#' @returns Prints the main results of the fitted APLMS model.
#' @examples
#' data(temperature)
#' datos = data.frame(temperature,time=1:length(temperature))
#' mod<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=datos,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' summary(mod)
#' @export
summary.aplms<-function(object, ...)
{
  summary_table <- generateSummaryTable(object)
  cat(" ---------------------------------------------------------------")
  cat("\n Additive partial linear models with symmetric errors \n")
  cat(" ---------------------------------------------------------------\n")
  cat(" Sample size: ", length(object$yhat))
  cat(" \n -------------------------- Model ---------------------------\n\n")
  print(object$this.call)
  cat(" \n ------------------- Parametric component -------------------\n\n")
  printCoefmat(summary_table, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat("\n ----------------- Non-parametric component ------------------ \n\n")
  printCoefmat(object$WALD_f, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(1,2))
  cat("\n --------------- Autoregressive and Scale parameter ---------------- \n\n")
  printCoefmat(object$summary_table_phirho, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat(" \n\n------ Penalized Log-likelihood and Information criterion------\n\n")
  cat(" Log-lik: ", round(object$loglike, digits = 2), "\n")
  cat(" AIC    : ", round(object$AIC, digits = 2), "\n")
  cat(" AICc   : ", round(object$AICC, digits = 2), "\n")
  cat(" BIC    : ", round(object$BIC, digits = 2), "\n")
  cat(" GCV    : ", round(object$GCV, digits = 2), "\n\n")
  cat(" --------------------------------------------------------------------\n")
}



# Generate the summary table of the parametric components.
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
generateSummaryTable <- function(model) {
  f0 <- model$f[[1]]
  est_coef <- as.vector(f0)
  ee <- sqrt(diag(model$VAR_F)[1:length(f0)])
  t_test <- est_coef / ee

  p_value <- sapply(t_test, FUN = function(x) {
    2 * pt(abs(x), model$rdf, lower.tail = FALSE)
  })

  summary_table <- cbind(
    est_coef,
    ee,
    t_test,
    p_value
  )

  terms_formula <- stats::terms(model$formula)
  var_names <- attr(terms_formula, "term.labels")

  rownames(summary_table) <- c("intercept", var_names)
  colnames(summary_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  return(summary_table)
}

# Generate the Walf F statistic and its table.
#
# @param npc_dimension cumulative sum of the number of rows of the function matrices
# @param dfk effective degrees of freedon
# @param npc vector with the names of non parametric components
# @param f estimated gamma parameters
# @param VAR_F covariance-variance matrix of the estimated gamma parameters
generateWaldF <- function(npc_dimension, dfk, npc, f, VAR_F) {
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

  WALD_f <- cbind(WALD_vec, dfk[-1], WALD_p_value)
  rownames(WALD_f) <- npc
  colnames(WALD_f) <- c("Wald", "df", "Pr(>.)")
  return(WALD_f)
}

# Generate the summary table of the rho parameter.
#
# @param p autoregressive order of the error
# @param par1 optimized log-likelihood function
# @param rho vector of estimated rho parameters
# @param rdf regressive degrees of freedom
generateSummaryTableRho <- function(p, par1, rho, rdf) {
  if (p > 0) {
    VAR_rho <- diag(solve(-par1$hessian))[2:(1 + p)]
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
      p_value_rho_t
    )

    colnames(summary_table_rho) <- c("rho", "ee", "Wald", "p-value_t")
    rownames(summary_table_rho) <- c(paste0("rho", as.character(1:p)))
  } else {
    summary_table_rho <- NULL
  }
  return(summary_table_rho)
}

# Generate the summary table of the rho parameter.
#
# @param p autoregressive order of the error
# @param par1 optimized log-likelihood function
# @param rdf regressive degrees of freedom
# @param rho vector of estimated rho parameters
# @param phi vector of estimated phi parameters
# @param nn number of observations.
# @param family probability distribution function.
generateSummartTablePhiRho <- function(p, par1, rdf, rho, phi, nn, family) {
  summary_table_rho <- generateSummaryTableRho(p, par1, rho, rdf)
  fg_t <- family$g3(args,
    df = family$df, r = family$r,
    s = family$s, alpha = family$alpha, mp = family$mp,
    epsi = family$epsi, sigmap = family$sigmap,
    k = family$k, nu = family$nu
  )
  fis_phi <- (nn / (4 * phi^2)) * (4 * fg_t - 1) # pag 52
  VAR_phi <- 1 / fis_phi

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
  return(summary_table_phirho)
}


