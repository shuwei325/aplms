#' @title Bayesian information criterion
#'
#' @description Print the BIC of the fitted APLMS model.
#' @param model APLMS object.
#'
#' @param ... other arguments
#' @examples
#' \dontrun{BIC(model)}
#' @method BIC aplms
#' @export

BIC.aplms<- function(model,...){
  if(!inherits(model, what="aplms", which = FALSE))
    stop("not a aplms object")
  model$BIC
}
