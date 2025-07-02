#' @title Akaike information criterion
#'
#' @description Print the AIC of the fitted APLMS model.
#' @param model APLMS object.
#'
#' @param ... other arguments
#' @examples
#' \dontrun{AIC(model)}
#' @method AIC aplms
#' @export

AIC.aplms<- function(model,...){
  if(!inherits(model, what="aplms", which = FALSE))
    stop("not a aplms object")
  model$AIC

}
