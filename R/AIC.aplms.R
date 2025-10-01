#' @title Akaike information criterion
#'
#' @description Print the AIC of the fitted APLMS model.
#' @param object APLMS object.
#' @param ... other arguments
#' @return a numeric value of the corresponding AIC.
#' @examples
#' \dontrun{AIC(object)}
#' @importFrom stats AIC
#' @method AIC aplms
#' @export

AIC.aplms<- function(object,...){
  if(!inherits(object, what="aplms", which = FALSE))
    stop("not a aplms object")
  object$AIC
}

