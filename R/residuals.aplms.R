#' Extract Residuals for APLMS fits
#'
#' @param object an object with the result of fitting additive partial linear models with symmetric errors.
#' \code{response} indicates response residuals, \code{pearson} is Pearson residuals, and \code{quant} is quantile residuals.
#' @param ... other arguments.
#' @return a dataframe with the residual, the pearson residual and the normal quantile of the standarized resiudals.
#' @keywords Additive partial linear models with symmetric errors
#' @keywords Residuals
#' @examples
#' \donttest{residuals(object)}
#' @importFrom stats residuals
#' @export
residuals.aplms <- function(object, ...) {
  if (!inherits(object, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }

  phi <- object$summary_table_phirho[1, 1]
  family_sym <- object$family

  xi_t <- family_sym$g4(1,
    df = family_sym$df,
    alpha = family_sym$alpha, mp = family_sym$mp, epsi = family_sym$epsi,
    sigmap = family_sym$sigmap, k = family_sym$k
  )
  residual_pearson <- object$residuals_y / sqrt(phi * xi_t)


  p_dist <- function(q, dist) {
    switch(dist,
      "Normal" = pnorm(q),
      "LogisI" = plogisI(q),
      "LogisII" = plogisII(q),
      "Student" = pt(q, df = family_sym$df),
      "Powerexp" = rmutil::ppowexp(q, m = 0, s = 1, f = 1 / (1 + family_sym$k)),
      "Gstudent" = pgstudent(q, s = family_sym$s, r = family_sym$r) # ,
      # 'Cauchy' =
      # 'Glogis' =
      # 'Cnormal' =
    )
  }

  res_stand <- object$residuals_y / sqrt(phi)
  residual_quant <- qnorm(p_dist(res_stand, family_sym$family))

  return(
    data.frame(
      res = as.numeric(object$residuals_y),
      res_pearson = as.numeric(residual_pearson),
      res_quant = as.numeric(residual_quant)
    )
  )
}


