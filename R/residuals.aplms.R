#' Extract Residuals for APLMS fits
#'
#' @param object an object with the result of fitting additive partial linear models with symmetric errors.
#' @param ... other arguments.
#' @return Returns a dataframe with the following columns
#' \describe{
#' \item{res}{the residual,}
#' \item{res_pearson}{the Pearson residual, and}
#' \item{res_quant}{the normal quantile of the standarized resiudals.}
#' }
#' @keywords Additive partial linear models with symmetric errors
#' @keywords Residuals
#' @examples
#' data(temperature)
#' temperature.df = data.frame(temperature,time=1:length(temperature))
#' model<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=temperature.df,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' residuals(model)
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
      "Gstudent" = pgstudent(q, s = family_sym$s, r = family_sym$r)
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


