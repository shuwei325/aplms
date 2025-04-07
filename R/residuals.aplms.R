# residuals.aplms(model,type=c("Quant","Pearson"),plot=NULL)


#' Extract Residuals for APLMS fits.
#'
#' @param mod an object with the result of fitting additive partial linear models with symmetric errors.
#' @param type a character string that indicates the type of residuals.
#' \code{response} indicates response residuals, \code{pearson} is Pearson residuals, and \code{quant} is quantile residuals.
#' @param plot a logical evaluating.
#' @param ... other arguments.
#' @keywords Additive partial linear models with symmetric errors
#' @keywords Residuals
#' @export residuals.aplms
#' @export
residuals.aplms <- function(mod, ...) {
  if (!inherits(mod, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }

  phi <- mod$summary_table_phirho[1, 1]
  family_sym <- mod$family

  xi_t <- family_sym$g4(1,
    df = family_sym$df,
    alpha = family_sym$alpha, mp = family_sym$mp, epsi = family_sym$epsi,
    sigmap = family_sym$sigmap, k = family_sym$k
  )
  residual_pearson <- mod$residuals_y / sqrt(phi * xi_t)


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
    )
  }

  res_stand <- mod$residuals_y / sqrt(phi)
  residual_quant <- qnorm(p_dist(res_stand, family_sym$family))


  # plots
  par(mfrow = c(3, 2))
  plot(residual_quant)
  hist(residual_quant)
  acf(residual_quant)
  pacf(residual_quant)
  qqnorm(residual_quant)
  qqline(residual_quant)

  return(
    data.frame(
      res = as.numeric(mod$residuals_y),
      res_pearson = as.numeric(residual_pearson),
      res_quant = as.numeric(residual_quant)
    )
  )
}
