#' Local influence plots of the object `aplms()`.
#'
#' Takes a fitted `aplms` object and outputs diagnostics of the sensitivity analysis by assessing the effects of perturbations in the model and/or data, on the parameter estimates. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
#' @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
#' @param C The cutoff criterion such that \eqn{C_i > \bar{C_i} + C*sd(C_i)} to detect influential observations.
#' @param labels label to especify each data point.
#' @export

influenceplot.aplms <- function(model,
                            perturbation = c("case-weight","dispersion","response","explanatory", "corAR"),
                            part = TRUE,
                            C = 4,
                            labels = NULL){

  perturbation <- match.arg(perturbation)
  if (is.null(labels)) labels <- 1:(nrow(model$data))

  output_list <- influence(model,
                    perturbation = perturbation,
                    part = part)

  if (perturbation %in% c("case-weight", "dispersion", "response")){

    print(paste0(perturbation," perturbation scheme"))

    if(part){
      par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
    } else {
      par(mfrow=c(1,1), oma = c(0, 0, 3, 0))
    }

    invisible(
      lapply(seq_along(output_list[[1]]),
           influential_plot1, output = output_list[[1]],
           labels = labels, C=C)
    )
    mtext(paste0(perturbation," perturbation scheme"), outer = TRUE, cex = 1.5)
  }


  if (perturbation %in% c("explanatory", "corAR")){

    print(paste0(perturbation," perturbation scheme"))

    old_par <- par(ask = TRUE)
    on.exit(par(old_par))

    if(part){
      par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
    } else {
      par(mfrow=c(1,1), oma = c(0, 0, 3, 0))
    }

    for( k in seq_along(output_list[[1]])) {
      invisible(
        lapply(seq_along(output_list[[1]][[k]]),
               influential_plot1, output = output_list[[1]][[k]],
               labels = labels, C=C)
        )
        mtext(paste0(perturbation," perturbation scheme - ", names(output_list[[1]])[k]),
              outer = TRUE, cex = 1.5)
    }
  }
}
