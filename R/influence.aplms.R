#' local influence analysis of the object `aplms()`
#'
#' Takes a fitted `aplms` object and outputs some diagnostic information about the fitting procedure and results. Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
#' @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
#' @param ... other arguments.
#' @return A list object containing the conformal normal curvature of the specified perturbations.
#' @examples
#' \dontrun{influence(model)}
#' @importFrom stats influence
#' @import psych MASS
#' @export
influence.aplms <- function(model,
                            perturbation = c("case-weight","dispersion","response","explanatory", "corAR"),
                            part = TRUE,...
){


  if (sapply(model$N_i,ncol)[1]==1 & any(perturbation == "explanatory")){
    stop("There is no explanatory variables.")
  }

  if ((nrow(model$summary_table_phirho)-1)== 0 & any(perturbation == "corAR")){
    stop("There is no autoregressive coeficients.")
  }


  list_output <- list()

  if (any(perturbation == "case-weight")){
    list_output <- perturbation_case1(model, list_output, "case-weight")
  }

  if (any(perturbation == "dispersion")){
    list_output <- perturbation_case1(model, list_output, "dispersion")
  }

  if (any(perturbation == "response")){
    list_output <- perturbation_case1(model, list_output, "response")
  }

  if (any(perturbation == "explanatory")){
    N_i <- model$N_i
    exp_index <- 1:(sapply(N_i,ncol)[1]-1)
    list_output[[which(perturbation=="explanatory")]] <- list()
    list_output <- perturbation_case2(model, list_output, "explanatory", exp_index)
    names(list_output[[which(perturbation=="explanatory")]]) <- row.names(model$f[[1]])[-1]
  }

  if (any(perturbation == "corAR")){
    p <- nrow(model$summary_table_phirho)-1
    k_index <- 1:p
    list_output[[which(perturbation=="corAR")]] <- list()
    list_output <- perturbation_case2(model, list_output, "corAR", k_index)
    names(list_output[[which(perturbation=="corAR")]]) <- paste0("rho",1:p)
  }

  names(list_output) <- perturbation
  return(list_output)
}

# Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response` perturbations are available.
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
# @param list_output an empty list object that will contain the conformal normal curvature of the specified perturbations.
# @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`.
# @return A list object containing the conformal normal curvature of the specified perturbations.
perturbation_case1 <- function(model, list_output, perturbation) {
  p <- nrow(model$summary_table_phirho)-1
  N_i <- model$N_i
  Lobs <- model$LL_obs
  coef_index <- c(sapply(N_i,dim)[2,],1,p)
  gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
  phi_index <- sum(coef_index)-p
  rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

  DELTA <- influence_DELTA(model = model, perturb_scheme = perturbation)
  C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)
  list_output[[which(perturbation==perturbation)]] <- list()
  list_output[[which(perturbation==perturbation)]]$C_i <- C_i
  if(part){
    # gamma
    Lobs.aux_gamma <- Gp(Lobs, index = gamma_index)
    C_i_gamma <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_gamma)
    # phi
    Lobs.aux_phi <- Gp(Lobs, index = phi_index)
    C_i_phi <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_phi)
    # rhos
    Lobs.aux_rhos <- Gp(Lobs, index = rho_index)
    C_i_rhos <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_rhos)
    list_output[[which(perturbation==perturbation)]]$C_i_gamma <- C_i_gamma
    list_output[[which(perturbation==perturbation)]]$C_i_phi <- C_i_phi
    list_output[[which(perturbation==perturbation)]]$C_i_rhos <- C_i_rhos
  }
  return(list_output)
}

# Returns the conformal normal curvature of the fitted `aplms` model object. The  `explanatory` and `corAR`. perturbations are available.
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
# @param list_output an empty list object that will contain the conformal normal curvature of the specified perturbations.
# @param perturbation A string vector specifying a perturbation scheme:  `explanatory` and `corAR`..
# @return A list object containing the conformal normal curvature of the specified perturbations.
perturbation_case2 <- function(model, list_output, perturbation, iteration) {
  p <- nrow(model$summary_table_phirho)-1
  N_i <- model$N_i
  Lobs <- model$LL_obs
  coef_index <- c(sapply(N_i,dim)[2,],1,p)
  gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
  phi_index <- sum(coef_index)-p
  rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

  for (i in iteration){
    if (perturbation == "explanatory") DELTA <- influence_DELTA(model = model, perturb_scheme = "explanatory", r = i)
    if (perturbation == "corAR") DELTA <- influence_DELTA(model = model, perturb_scheme = "corAR", k = i)
    C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)
    list_output[[which(perturbation==perturbation)]][[i]] <- list()
    list_output[[which(perturbation==perturbation)]][[i]]$C_i <- C_i
    if(part){
      # gamma
      Lobs.aux_gamma <- Gp(Lobs, index = gamma_index)
      C_i_gamma <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_gamma)
      # phi
      Lobs.aux_phi <- Gp(Lobs, index = phi_index)
      C_i_phi <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_phi)
      # rhos
      Lobs.aux_rhos <- Gp(Lobs, index = rho_index)
      C_i_rhos <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_rhos)
      list_output[[which(perturbation==perturbation)]][[i]]$C_i_gamma <- C_i_gamma
      list_output[[which(perturbation==perturbation)]][[i]]$C_i_phi <- C_i_phi
      list_output[[which(perturbation==perturbation)]][[i]]$C_i_rhos <- C_i_rhos
    }
  }
  return(list_output)
}