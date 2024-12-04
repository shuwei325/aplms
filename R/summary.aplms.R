#' Print method for "foo" class
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param ... Other arguments passed to or from other methods
#' @rdname summary
#' @method summary aplms
#' @export summary.aplms
#' @export
summary.aplms<-function(model, ...)
{
  cat(" ---------------------------------------------------------------")
  cat("\n Additive partial linear models with symmetric errors \n")
  cat(" ---------------------------------------------------------------\n")
  cat(" Sample size: ", length(model$yhat))
  cat(" \n -------------------------- Model ---------------------------\n\n")
  print(model$this.call)
  cat(" \n ------------------- Parametric component -------------------\n\n")
  printCoefmat(model$summary_table, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat("\n ----------------- Non-parametric component ------------------ \n\n")
  printCoefmat(model$WALD_f, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(1,2))
  cat("\n --------------- Autoregressive and Scale parameter ---------------- \n\n")
  printCoefmat(model$summary_table_phirho, P.values = TRUE, has.Pvalue = TRUE, digits = 5,
               signif.legend = FALSE, tst.ind = c(2, 3))
  cat(" \n\n------ Penalized Log-likelihood and Information criterion------\n\n")
  cat(" Log-lik: ", round(model$loglike, digits = 2), "\n")
  cat(" AIC    : ", round(model$AIC, digits = 2), "\n")
  cat(" AICc   : ", round(model$AICC, digits = 2), "\n")
  cat(" BIC    : ", round(model$BIC, digits = 2), "\n")
  cat(" GCV    : ", round(model$GCV, digits = 2), "\n\n")
  cat(" --------------------------------------------------------------------\n")
}
