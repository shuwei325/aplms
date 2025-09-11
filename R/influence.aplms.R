#' local influence analysis of the object `aplms()`.
#'
#' Takes a fitted `aplms` object and outputs some diagnostic information about the fitting procedure and results. Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
#' @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
#' @examples
#' \dontrun{influence(model)}
#' @import psych MASS
#' @method influence aplms
#' @export

influence.aplms <- function(model,
                            perturbation = c("case-weight","dispersion","response","explanatory", "corAR"),
                            part = TRUE
){


  if (sapply(model$N_i,ncol)[1]==1 & any(perturbation == "explanatory")){
    stop("There is no explanatory variables.")
  }

  if ((nrow(model$summary_table_phirho)-1)== 0 & any(perturbation == "corAR")){
    stop("There is no autoregressive coeficients.")
  }


  list_output <- list()

   p <- nrow(model$summary_table_phirho)-1
  # nn <- nrow(model$data)
  # ONE <- cbind(rep(1,nn))
  # A <- matrix_A(model$summary_table_phirho[2:(p+1),1],nn)
  # B <- BB(p=p,nn)
  #
  # phi <- model$summary_table_phirho[1,1]
  N_i <- model$N_i
  Lobs <- model$LL_obs

  coef_index <- c(sapply(N_i,dim)[2,],1,p)
  gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
  phi_index <- sum(coef_index)-p
  rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

  if (any(perturbation == "case-weight")){

    DELTA <- influence_DELTA(model = model, perturb_scheme = "case-weight")

    C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)

    list_output[[which(perturbation=="case-weight")]] <- list()
    list_output[[which(perturbation=="case-weight")]]$C_i <- C_i

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

      list_output[[which(perturbation=="case-weight")]]$C_i_gamma <- C_i_gamma
      list_output[[which(perturbation=="case-weight")]]$C_i_phi <- C_i_phi
      list_output[[which(perturbation=="case-weight")]]$C_i_rhos <- C_i_rhos

    }
  }

  if (any(perturbation == "dispersion")){

    DELTA <- influence_DELTA(model = model, perturb_scheme = "dispersion")
    C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)

    list_output[[which(perturbation=="dispersion")]] <- list()
    list_output[[which(perturbation=="dispersion")]]$C_i <- C_i

    if(part){

      # gamma
      Lobs.aux_gamma <- Gp(Lobs, index = gamma_index)
      C_i_gamma <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_gamma)

      # phi
      Lobs.aux_phi <- Gp(Lobs, index = phi_index)
      C_i_phi <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_phi)

      # rhos
      Lobs.aux_rhos <- Gp(Lobs, index = rho_index)
      C_i_rhos <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_rhos)

      list_output[[which(perturbation=="dispersion")]]$C_i_gamma <- C_i_gamma
      list_output[[which(perturbation=="dispersion")]]$C_i_phi <- C_i_phi
      list_output[[which(perturbation=="dispersion")]]$C_i_rhos <- C_i_rhos
    }
  }

  if (any(perturbation == "response")){

    DELTA <- influence_DELTA(model = model, perturb_scheme = "response")

    C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)

    list_output[[which(perturbation=="response")]] <- list()
    list_output[[which(perturbation=="response")]]$C_i <- C_i

    if(part){

      # gamma
      Lobs.aux_gamma <- Gp(Lobs, index = gamma_index)
      C_i_gamma <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_gamma)

      # phi
      Lobs.aux_phi <- Gp(Lobs, index = phi_index)
      C_i_phi <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_phi)

      # rhos
      Lobs.aux_rhos <- Gp(Lobs, index = rho_index)
      C_i_rhos <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = Lobs.aux_rhos)

      list_output[[which(perturbation=="response")]]$C_i_gamma <- C_i_gamma
      list_output[[which(perturbation=="response")]]$C_i_phi <- C_i_phi
      list_output[[which(perturbation=="response")]]$C_i_rhos <- C_i_rhos
    }

  }

  if (any(perturbation == "explanatory")){

      exp_index <- 1:(sapply(N_i,ncol)[1]-1)
      list_output[[which(perturbation=="explanatory")]] <- list()

      for (r in exp_index){

        DELTA <- influence_DELTA(model = model, perturb_scheme = "explanatory", r = r)

        C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)

        list_output[[which(perturbation=="explanatory")]][[r]] <- list()
        list_output[[which(perturbation=="explanatory")]][[r]]$C_i <- C_i

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

          list_output[[which(perturbation=="explanatory")]][[r]]$C_i_gamma <- C_i_gamma
          list_output[[which(perturbation=="explanatory")]][[r]]$C_i_phi <- C_i_phi
          list_output[[which(perturbation=="explanatory")]][[r]]$C_i_rhos <- C_i_rhos
        }
      }
    names(list_output[[which(perturbation=="explanatory")]]) <- row.names(model$f[[1]])[-1]
  }

  if (any(perturbation == "corAR")){

      k_index <- 1:p

      list_output[[which(perturbation=="corAR")]] <- list()

      for (k in k_index){

        DELTA <- influence_DELTA(model = model, perturb_scheme = "corAR", k = k)

        C_i <- conf_normal_curvature(DELTA, Lobs, Lobs.aux = 0)

        list_output[[which(perturbation=="corAR")]][[k]] <- list()
        list_output[[which(perturbation=="corAR")]][[k]]$C_i <- C_i

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

          list_output[[which(perturbation=="corAR")]][[k]]$C_i_gamma <- C_i_gamma
          list_output[[which(perturbation=="corAR")]][[k]]$C_i_phi <- C_i_phi
          list_output[[which(perturbation=="corAR")]][[k]]$C_i_rhos <- C_i_rhos
        }
      }


    names(list_output[[which(perturbation=="corAR")]]) <- paste0("rho",1:p)

  }

  names(list_output) <- perturbation

  return(list_output)

}




