#' Diagnostic Plots for additive partial linear models with symmetric errors
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param which an optional numeric value with the number of only plot that must be returned.
#' @param labels a optional string vector specifying a labels plots.
#' @param iden a logical value used to identify observations. If \code{TRUE} the observations are identified by user in the graphic window.
#' @param ... graphics parameters to be passed to the plotting routines.
#' @return Return an interactive menu with eleven options to make plots. This menu contains the follows graphics:
#' 1: Response residuals against fited values.
#' 2: Response residuals against time index.
#' 3: Histogram of Response residuals.
#' 4: Autocorrelation function of response residuals.
#' 5: Partial autocorrelation function of response residuals.
#' 6: Conditional quantile residuals against fited values.
#' 7: Conditional quantile residuals against time index.
#' 8: Histogram of conditional quantile residuals.
#' 9: Autocorrelation function of conditional quantile residual.
#' 10: Partial autocorrelation function of conditional quantile residuals.
#' 11: QQ-plot of conditional quantile residuals.
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
#' aplms.diag.plot(model, which = 1)
#' @keywords Additive partial linear models with symmetric errors
#' @keywords Residuals
#' @import utils graphics
#' @export
aplms.diag.plot <- function(model, which,
                            labels = NULL,
                            iden = FALSE, ...) {
  if (!inherits(model, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }

  res_data <- residuals(model)

  choices <- c(
               "Response residuals against fited values",
               "Response residuals against time index",
               "Histogram of Response residuals",
               "Autocorrelation function of response residuals",
               "Partial autocorrelation function of response residuals",
               "Conditional quantile residuals against fited values",
               "Conditional quantile residuals against time index",
               "Histogram of conditional quantile residuals",
               "Autocorrelation function of conditional quantile residuals",
               "Partial autocorrelation function of conditional quantile residuals",
               "QQ-plot of conditional quantile residuals"
               )
  n_choice <- seq_along(choices)


  tmenu <- paste("plot:", choices)
  if (missing(which))
    pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
  else if (!match(which, n_choice, nomatch = F))
    stop("choice not valid")
  else pick <- which
  if (pick == 0)
    stop(" no graph required ! ")

  repeat {
    switch(pick,
           `1` = {
             y1 <- res_data$res
             x1 <- model$yhat
            plot(x1, y1, xlab = "Fitted values", ylab = "Response residual",
                  ...)
            grid()
            xx <- list(x1)
            yy <- list(y1)
           },
           `2` = {
             y2 <- res_data$res
             x2 <- seq_along(y2)
             plot(x2, y2, xlab = "Index", ylab = "Response residual",
                  ...)
             grid()
             xx <- list(x2)
             yy <- list(y2)
           },
           `3` = {
             y3 <- res_data$res
             hist(y3, main="Response residuals")
             grid()
           },
           `4` = {
             y4 <- res_data$res
             acf(y4, main="Response residuals")
           },
           `5` = {
             y5 <- res_data$res
             pacf(y5, main="Response residuals")
           },
           `6` = {
             y6 <- res_data$res_quant
             x6 <- model$yhat
             plot(x6, y6, xlab = "Fitted values", ylab = "Conditional quantile residual",
                  ...)
             grid()
             xx <- list(x6)
             yy <- list(y6)
           },
           `7` = {
             y7 <- res_data$res_quant
             x7 <- seq_along(y7)
             plot(x7, y7, xlab = "Index", ylab = "Conditional quantile residual",
                  ...)
             grid()
             xx <- list(x7)
             yy <- list(y7)
           },
           `8` = {
             y8 <- res_data$res_quant
             hist(y8, main="Conditional quantile residual")
             grid()
           },
           `9` = {
             y9 <- res_data$res_quant
             acf(y9, main="Conditional quantile residual")
           },
           `10` = {
             y10<- res_data$res_quant
             pacf(y10, main="Conditional quantile residual")
           },
           `11` = {
             y11 <- res_data$res_quant
             x11 <- qnorm(ppoints(length(y11)))[rank(y11)]
             .lim <- c(min(x11, y11), max(x11, y11))
             plot(x11, y11, xlab = paste("Quantiles of standard normal"),
                  ylab = "Ordered standardized residual", xlim = .lim,
                  ylim = .lim, ...)
             grid()
             abline(0, 1, lty = 2)
             xx <- list(x11)
             yy <- list(y11)
           })

    if (pick == 1 || pick == 2 || pick == 6 || pick == 7 || pick == 11 ) {
      if (is.null(labels))
        labels <- seq_along(res_data$res)
      yes <- iden
      while (yes) {
        cat("****************************************************\n")
        cat("Interactive Identification\n")
        cat("left button = Identify, center button = Exit\n")
        identify(xx[[1]], yy[[1]], labels, ...)
        yes <- F
      }
    }
    if (missing(which))
      pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
    if ((pick == 0) || !missing(which)) {
      invisible(close.screen(all.screens = T))
      break
    }
  }
}
