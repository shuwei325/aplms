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
#' data(temperature)
#' temperature.df = data.frame(temperature,time=1:length(temperature))
#' model<-aplms::aplms(temperature ~ 1,
#'                    npc=c("time"), basis=c("cr"),Knot=c(60),
#'                    data=temperature.df,family=Powerexp(k=0.3),p=1,
#'                    control = list(tol = 0.001,
#'                                   algorithm1 = c("P-GAM"),
#'                                   algorithm2 = c("BFGS"),
#'                                   Maxiter1 = 20,
#'                                   Maxiter2 = 25),
#'                    lam=c(10))
#' influence(model, perturbation = c("case-weight"))
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
    list_output[[which(perturbation=="case-weight")]] <- perturbation_case1(model, "case-weight", part)
  }

  if (any(perturbation == "dispersion")){
    list_output[[which(perturbation=="dispersion")]] <- perturbation_case1(model, "dispersion", part)
  }

  if (any(perturbation == "response")){
    list_output[[which(perturbation=="response")]] <- perturbation_case1(model, "response", part)
  }

  if (any(perturbation == "explanatory")){
    N_i <- model$N_i
    exp_index <- 1:(sapply(N_i,ncol)[1]-1)
    list_output[[which(perturbation=="explanatory")]] <- perturbation_case2(model, "explanatory", exp_index, part)
    names(list_output[[which(perturbation=="explanatory")]]) <- attr(stats::terms(model$formula), "term.labels")
  }

  if (any(perturbation == "corAR")){
    p <- nrow(model$summary_table_phirho)-1
    k_index <- 1:p
    list_output[[which(perturbation=="corAR")]] <- perturbation_case2(model, "corAR", k_index, part)
    names(list_output[[which(perturbation=="corAR")]]) <- paste0("rho",1:p)
  }

  names(list_output) <- perturbation
  return(list_output)
}

# Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response` perturbations are available.
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
# @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`.
# @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
# @return A list object containing the conformal normal curvature of the specified perturbations.
perturbation_case1 <- function(model, perturbation, part) {
  p <- nrow(model$summary_table_phirho)-1
  N_i <- model$N_i
  Lobs <- model$LL_obs
  coef_index <- c(sapply(N_i,dim)[2,],1,p)
  gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
  phi_index <- sum(coef_index)-p
  rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

  DELTA <- influence_DELTA(model = model, perturb_scheme = perturbation)
  C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)
  list_output <- list()
  list_output$C_i <- C_i
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
    list_output$C_i_gamma <- C_i_gamma
    list_output$C_i_phi <- C_i_phi
    list_output$C_i_rhos <- C_i_rhos
  }
  return(list_output)
}

# Returns the conformal normal curvature of the fitted `aplms` model object. The  `explanatory` and `corAR`. perturbations are available.
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
# @param perturbation A string vector specifying a perturbation scheme:  `explanatory` and `corAR`.
# @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
# @return A list object containing the conformal normal curvature of the specified perturbations.
perturbation_case2 <- function(model, perturbation, iteration, part) {
  p <- nrow(model$summary_table_phirho)-1
  N_i <- model$N_i
  Lobs <- model$LL_obs
  coef_index <- c(sapply(N_i,dim)[2,],1,p)
  gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
  phi_index <- sum(coef_index)-p
  rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

  list_output <- list()
  for (i in iteration){
    if (perturbation == "explanatory") DELTA <- influence_DELTA(model = model, perturb_scheme = "explanatory", r = i)
    if (perturbation == "corAR") DELTA <- influence_DELTA(model = model, perturb_scheme = "corAR", k = i)
    C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)
    list_output[[i]] <- list()
    list_output[[i]]$C_i <- C_i
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
      list_output[[i]]$C_i_gamma <- C_i_gamma
      list_output[[i]]$C_i_phi <- C_i_phi
      list_output[[i]]$C_i_rhos <- C_i_rhos
    }
  }
  return(list_output)
}
