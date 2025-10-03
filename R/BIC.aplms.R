#' @title Bayesian information criterion
#'
#' @description Print the BIC of the fitted APLMS model.
#' @param object APLMS object.
#' @param ... other arguments
#' @return a numeric value of the corresponding BIC.
#' @examples
#' \dontrun{BIC(object)}
#' @importFrom stats BIC
#' @export
BIC.aplms<- function(object,...){
  if(!inherits(object, what="aplms", which = FALSE))
    stop("not a aplms object")
  object$BIC
}

