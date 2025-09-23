#' @title Bayesian information criterion
#'
#' @description Print the BIC of the fitted APLMS model.
#' @param model APLMS object.
#' @param ... other arguments
#' @return a numeric value of the corresponding BIC.
#' @examples
#' \dontrun{BIC(model)}
#' @method BIC aplms
#' @export
BIC.aplms<- function(model,...){
  if(!inherits(model, what="aplms", which = FALSE))
    stop("not a aplms object")
  model$BIC
}

#' @rdname BIC.aplms
#' @export
BIC <- function(model, ...) {
  UseMethod("BIC")
}
