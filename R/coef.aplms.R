#' Extract the coefficients of the fitted APLMS model
#'
#' Extract the coefficients of the fitted APLMS model.
#' @param object APLMS object.
#' @param ... other arguments
#' @return A list of vectors of the corresponding estimated parameters.
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
#' coef(model)
#' @importFrom stats coef
#' @method coef aplms
#' @export
coef.aplms<- function(object,...){
  if(!inherits(object, what="aplms", which = FALSE))
    stop("not a aplms object")
  list(
    gamma = object$summary_table[,1],
    phirho = object$WALD_f[,1])

}


