#' Local influence plots of the object `aplms()`
#'
#' Takes a fitted `aplms` object and outputs diagnostics of the sensitivity analysis by assessing the effects of perturbations in the model and/or data, on the parameter estimates. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
#' @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
#' @param C The cutoff criterion such that \eqn{C_i > \bar{C_i} + C*sd(C_i)} to detect influential observations.
#' @param labels label to especify each data point.
#' @return The conformal normal curvature of the specified perturbations is plotted.
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
#' influenceplot.aplms(model, perturbation = c("case-weight"))
#' @importFrom grid gpar
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 geom_point geom_hline geom_text
#' @export

influenceplot.aplms <- function(model,
                            perturbation = c("case-weight","dispersion","response","explanatory", "corAR"),
                            part = TRUE,
                            C = 4,
                            labels = NULL){

  perturbation <- match.arg(perturbation)

  if (!(is.numeric(C) & C>0)) stop("C should be a positive value.")

  if (!is.null(labels)) {
    if (length(labels) != 2) {
      stop("labels length should be equal to number of observations.")
    }
  } else {
    labels <- seq_along(model$yhat)
  }

  output_list <- influence(model,
                    perturbation = perturbation,
                    part = part)


  if (perturbation %in% c("case-weight", "dispersion", "response")){

    print(paste0(perturbation," perturbation scheme"))

    if(part){
      nrow <- 2
      ncol <- 2
    } else {
      nrow <- 1
      ncol <- 1
    }

    invisible(
      plots <- lapply(seq_along(output_list[[1]]),
           influential_plot1, output = output_list[[1]],
           labels = labels, C=C)
    )
    grid.arrange(grobs = plots, nrow = nrow, ncol = ncol,
             top = paste0(perturbation," perturbation scheme"), gp = gpar(fontsize = 24, fontface = "bold"))
  }


  if (perturbation %in% c("explanatory", "corAR")){

    print(paste0(perturbation," perturbation scheme"))

    if(part){
      nrow <- 2
      ncol <- 2
    } else {
      nrow <- 1
      ncol <- 1
    }

    for( k in seq_along(output_list[[1]])) {
      invisible(
        plots <- lapply(seq_along(output_list[[1]][[k]]),
               influential_plot1, output = output_list[[1]][[k]],
               labels = labels, C=C)
        )
        grid.arrange(grobs = plots, nrow = nrow, ncol = ncol,
              top = paste0(perturbation," perturbation scheme - ", names(output_list[[1]])[k]),
              gp = gpar(fontsize = 24, fontface = "bold"))
    }
  }
}
