#' @title Bayesian information criterion
#'
#' @description Print the BIC of the fitted APLMS model.
#' @param object APLMS object.
#' @param ... other arguments
#' @return a numeric value of the corresponding BIC.
#' @examples
#' \donttest{
#' data(temperature)
#' datos = data.frame(temperature,time=1:length(temperature))
#' mod<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=datos,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' BIC(mod)
#' }
#' @importFrom stats BIC
#' @export
BIC.aplms<- function(object,...){
  if(!inherits(object, what="aplms", which = FALSE))
    stop("not a aplms object")
  object$BIC
}

