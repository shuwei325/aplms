#' Extract the coefficients of the fitted APLMS model
#'
#' Extract the coefficients of the fitted APLMS model.
#' @param model APLMS object.
#' @param ... other arguments
#' @return A list vector of the corresponding estimated parameters.
#' @examples
#' \dontrun{coef(model)}
#' @method coef aplms
#' @export
coef.aplms<- function(model,...){
  if(!inherits(model, what="aplms", which = FALSE))
    stop("not a aplms object")
  list(
    gamma = model$summary_table[,1],
    phirho = model$WALD_f[,1])

}

#' @rdname coef.aplms
#' @export
coef <- function(model, ...) {
  UseMethod("coef")
}
