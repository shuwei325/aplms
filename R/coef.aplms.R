#' Extract the coefficients of the fitted APLMS model
#'
#' Extract the coefficients of the fitted APLMS model.
#' @param object APLMS object.
#' @param ... other arguments
#' @return A list vector of the corresponding estimated parameters.
#' @examples
#' \dontrun{coef(model)}
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


