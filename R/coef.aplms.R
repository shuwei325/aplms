#' Extract the coefficients of the fitted APLMS model
#'
#' Extract the coefficients of the fitted APLMS model.
#' @param model APLMS object.
#'
#' @param ... other arguments
#' @examples
#' a <- 1
#' @method coef aplms
#' @export

coef.aplms<- function(model,...){
  if(!inherits(mod, what="aplms", which = FALSE))
    stop("not a aplms object")
  list(
    gamma = model$summary_table[,1],
    phirho = model$WALD_f[,1])

}
