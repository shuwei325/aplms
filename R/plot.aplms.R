#' Default APLMS plotting.
#'
#' Compute and plot the estimated mean and confidence intervals of the non-parametric component of a `APLMS` object fited by `aplms()`.
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param level Confidence level.
#' @param len The desired length of the sequence of covariates to compute the non parametric component functions.
#' @param plot a logical value to return plots. Default value is \code{TRUE}.
#' @param ... other arguments.
#' @return Return a list of all non parametric component functions with their confidence intervals.
#'   interactive menu with eleven options to make plots.
#' @keywords Additive partial linear models with symmetric errors
#' @keywords Residuals
#' @examples
#' \dontrun{
#' data(temperature)
#' datos = data.frame(temperature,time=1:length(temperature))
#' mod1<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=datos,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' plot(mod1)
#' }
#' @method plot aplms
#' @export
plot.aplms <- function(model, level = 0.95, len = 100, plot = TRUE , ...) {
  if (!inherits(model, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }

  npc_predict_list <- list()

  for (k in seq_along(model$npc)) {
    cov_name <- model$npc[(k)]
    basis <- model$basis[(k)]
    f_cov <- model$f[[(k+1)]]
    nm <- model$Knot[(k)]

    data_aux<-data.frame(cov=seq(min(model$data[cov_name]),max(model$data[cov_name]),length = len))
    ZZ_aux<-smoothCon(s(cov,bs=basis,k=nm),data=data_aux,knots=NULL,absorb.cons=T)
    N1_aux<-ZZ_aux[[1]]$X

    fmean <- N1_aux %*% f_cov
    VAR_F <- model$VAR_F
    dim_coef <- sapply(model$N_i,dim)[2,]
    index_dim_coef <- cumsum(dim_coef)

    se <- sqrt(
      diag(
        N1_aux%*%tcrossprod(VAR_F[
          (1+index_dim_coef[(k)]):(index_dim_coef[(k+1)]),
          (1+index_dim_coef[(k)]):(index_dim_coef[(k+1)])],
          N1_aux)
      )
    )

    data_aux$fmean = fmean
    data_aux$fmean_ls = fmean+qnorm((1-level)/2/nm)*se
    data_aux$fmean_li = fmean-qnorm((1-level)/2/nm)*se

    npc_predict_list[[k]] <- data_aux
  }

  if (plot) {

    if (length(model$npc)==1) {
      split.screen(c(1, 1))
    } else {
      split.screen(c(ceiling(length(model$npc)/2), 2))
    }

    for (k in seq_along(model$npc)) {
      screen(k)
      plot(npc_predict_list[[k]]$cov,
           npc_predict_list[[k]]$fmean,
           type = "l",
           ylim = range(c(npc_predict_list[[k]]$fmean_li, npc_predict_list[[k]]$fmean_ls)),
           xlab = model$npc[k], ylab = "effect", main = model$npc[k])
      points(npc_predict_list[[k]]$cov,
             npc_predict_list[[k]]$fmean_li, type = "l", lty = 2)
      points(npc_predict_list[[k]]$cov,
             npc_predict_list[[k]]$fmean_ls, type = "l", lty = 2)
    }
    close.screen(all.screens = TRUE)
  }

return(npc_predict_list)

}

#' @rdname plot.aplms
#' @export
plot.aplms <- function(model, level, len, plot, ...) {
  UseMethod("plot")
}
