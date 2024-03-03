#' Print method for "foo" class
#'
#' @param object an object with the result of fitting additive partial linear models with symmetric errors.
#' @param ... Other arguments passed to or from other methods
#'
#' @export
summary  <- function(object, ...) {
  UseMethod("summary")
}

#' @rdname summary
#' @export
summary.aplms<-function(object, ...)
{
  cat(" ---------------------------------------------------------------")
  cat("\n Additive partial linear models with symmetric errors \n")
  cat(" ---------------------------------------------------------------\n")
  cat(" Sample size: ", length(object$yhat))
  cat(" \n -------------------------- Model ---------------------------\n\n")
  cat(" \n ------------------- Parametric component -------------------\n\n")
  printCoefmat(object$summary_table, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat("\n ----------------- Non-parametric component ------------------ \n\n")
  printCoefmat(object$WALD_f, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(1,2))
  cat("\n --------------- Autoregressive and Scale parameter ---------------- \n\n")
  printCoefmat(object$summary_table_rhophi, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat(" \n\n------ Penalized Log-likelihood and Information criterion------\n\n")
  cat(" Log-lik: ", round(object$loglike, digits = 2), "\n")
  cat(" AIC    : ", round(object$AIC, digits = 2), "\n")
  cat(" AICc   : ", round(object$AICC, digits = 2), "\n")
  cat(" BIC    : ", round(object$BIC, digits = 2), "\n")
  cat(" GCV    : ", round(object$GCV, digits = 2), "\n\n")
  cat(" --------------------------------------------------------------------\n")
}
