#' Default APLMS plotting
#'
#' Compute and plot the estimated mean and confidence intervals of the non-parametric component of a `APLMS` object fited by `aplms()`.
#' @param x an object with the result of fitting additive partial linear models with symmetric errors.
#' @param len The desired length of the sequence of covariates to compute the non parametric component functions.
#' @param plot a logical value to return plots. Default value is \code{TRUE}.
#' @param level Confidence level.
#' @param ... other arguments.
#' @return Return a list of all non parametric component functions with their confidence intervals. If \code{plot=TRUE}, the estimated nonparametric component functions are plotted.
#' @keywords Additive partial linear xs with symmetric errors
#' @keywords Residuals
#' @importFrom graphics plot
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
#' plot(model)
#' @export
plot.aplms <- function(x, len = 100, plot = TRUE , level = 0.95, ...) {
  if (!inherits(x, what = "aplms", which = FALSE)) {
    stop("not a aplms object")
  }

  npc_predict_list <- list()

  for (k in seq_along(x$npc)) {
    cov_name <- x$npc[(k)]
    basis <- x$basis[(k)]
    f_cov <- x$f[[(k+1)]]
    nm <- x$Knot[(k)]

    data_aux<-data.frame(cov=seq(min(x$data[cov_name]),max(x$data[cov_name]),length = len))
    ZZ_aux<-smoothCon(s(cov,bs=basis,k=nm),data=data_aux,knots=NULL,absorb.cons=T)
    N1_aux<-ZZ_aux[[1]]$X

    fmean <- N1_aux %*% f_cov
    VAR_F <- x$VAR_F
    dim_coef <- sapply(x$N_i,dim)[2,]
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

    if (length(x$npc)==1) {
      split.screen(c(1, 1))
    } else {
      split.screen(c(ceiling(length(x$npc)/2), 2))
    }

    for (k in seq_along(x$npc)) {
      screen(k)
      plot(npc_predict_list[[k]]$cov,
           npc_predict_list[[k]]$fmean,
           type = "l",
           ylim = range(c(npc_predict_list[[k]]$fmean_li, npc_predict_list[[k]]$fmean_ls)),
           xlab = x$npc[k], ylab = parse(text = paste0("f[", k,"](.)")), main = x$npc[k])
      points(npc_predict_list[[k]]$cov,
             npc_predict_list[[k]]$fmean_li, type = "l", lty = 2)
      points(npc_predict_list[[k]]$cov,
             npc_predict_list[[k]]$fmean_ls, type = "l", lty = 2)
      grid()
    }
    close.screen(all.screens = TRUE)
  }

return(npc_predict_list)
}


