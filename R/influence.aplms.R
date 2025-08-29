#' local influence of the object `aplms()`.
#'
#' Takes a fitted `aplms` object and outputs some diagnostic information about the fitting procedure and results. Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
#' @param part A logical value to indicate whether the influential analysis is performed for \eqn{\gamma}, \eqn{\phi} and \eqn{\rho}.
#' @param C The cutoff criterion such that \eqn{C_i > \bar{C_i} + C*sd(C_i)} to detect influential observations.
#' @param plot A logical value used to plot.
#' @import psych MASS
#' @export

influence.aplms <- function(model,
                perturbation = c("case-weight","dispersion","response","explanatory"), #,, "corAR"),
                part = TRUE,
                C = 4,
                #ident = NULL,
                plot = FALSE#,
                # labels = NULL,
                # iden = F
                ){

list_output <- list()

p <- nrow(model$summary_table_phirho)-1
nn <- nrow(model$data)
ONE <- cbind(rep(1,nn))
A <- matrix_A(model$summary_table_phirho[2:(p+1),1],nn)
B <- BB(p=p,nn)

phi <- model$summary_table_phirho[1,1]
N_i <- model$N_i
Lobs <- model$LL_obs

coef_index <- c(sapply(N_i,dim)[2,],1,p)
gamma_index <- seq_len(sum(sapply(N_i,dim)[2,]))
phi_index <- sum(coef_index)-p
rho_index <- ((sum(coef_index)-p+1):sum(coef_index))

if (any(perturbation == "case-weight")){

    #calculo de las Deltas
    Delta_i <- lapply(model$N_i,
                      FUN = function(x){
                          t(t(A%*%x) %*% model$Dv %*% diag(c(A %*% model$residuals_mu))/phi)
                        })
    Delta_phi<- -ONE/(2*phi) + model$Dm%*%ONE/(2*phi)
    Delta_rho <- lapply(B,
                        FUN = function(x){t(t(x %*% model$residuals_mu) %*%
                                              model$Dv %*% diag(c(A %*% model$residuals_mu))/phi)})
    DELTA <- t(cbind(do.call(cbind,Delta_i),
                     Delta_phi,
                     do.call(cbind,Delta_rho)))

    C_i <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = 0)
    C_i_iden <- influential_index(C_i, C = 4)

    list_output[[which(perturbation=="case-weight")]] <- list()
    list_output[[which(perturbation=="case-weight")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    # gamma

    B_gamma <- Lobs
    B_gamma[gamma_index,] <- 0
    B_gamma[,gamma_index] <- 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)

    C_i_gamma <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_gamma)
    C_i_gamma_iden <- influential_index(C_i_gamma, C= 4)

    # phi

    B_phi <- Lobs
    B_phi[phi_index,] <- 0
    B_phi[,phi_index] <- 0
    Lobs.aux_phi = MASS::ginv(B_phi)

    C_i_phi <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_phi)
    C_i_phi_iden <- influential_index(C_i_phi, C= 4)

    # rhos

    B_rhos <- Lobs
    B_rhos[rho_index,]<- 0
    B_rhos[,rho_index] <- 0
    Lobs.aux_rhos = MASS::ginv(B_rhos)

    C_i_rhos <- conf_normal_curvature(DELTA,Lobs, Lobs.aux = Lobs.aux_rhos)
    C_i_rhos_iden <- influential_index(C_i_rhos, C= 4)

    list_output[[which(perturbation=="case-weight")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                 C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                 C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }
}

if (any(perturbation == "dispersion")){

  #calculo de las Deltas

  Delta_i <- lapply(N_i,
                    FUN = function(x){t(t(A%*%x) %*% (model$Dv-2*model$Dd) %*%
                                          diag(c(A %*% model$residuals_mu))/phi)})
  Delta_phi<- 1/(2*phi)*t(model$delta) %*% ( model$Dv - 2* model$Dd)
  Delta_rho <- lapply(B,
                      FUN = function(x){t( t(x%*%model$residuals_mu)%*%
                                             (2*model$Dd-model$Dv)%*%
                                             diag(c(A%*%model$residuals_mu)) /phi)})
  DELTA_pdispersao <- t(cbind(do.call(cbind,Delta_i),
                              t(Delta_phi),
                              do.call(cbind,Delta_rho)))

  C_i <- conf_normal_curvature(DELTA_pdispersao, Lobs, Lobs.aux = 0)
  C_i_iden <- influential_index(C_i, C = 4)


  list_output[[which(perturbation=="dispersion")]] <- list()
  list_output[[which(perturbation=="dispersion")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    # gamma

    B_gamma <- Lobs
    B_gamma[gamma_index,] <- 0
    B_gamma[,gamma_index] <- 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)

    C_i_gamma <- conf_normal_curvature(DELTA_pdispersao, Lobs, Lobs.aux = Lobs.aux_gamma)
    C_i_gamma_iden <- influential_index(C_i_gamma, C = 4)

    # phi

    B_phi <- Lobs
    B_phi[phi_index,] <- 0
    B_phi[,phi_index] <- 0
    Lobs.aux_phi = MASS::ginv(B_phi)

    C_i_phi <- conf_normal_curvature(DELTA_pdispersao, Lobs, Lobs.aux = Lobs.aux_phi)
    C_i_phi_iden <- influential_index(C_i_phi, C = 4)

    # rhos

    B_rhos <- Lobs
    B_rhos[rho_index,]<- 0
    B_rhos[,rho_index]<- 0
    Lobs.aux_rhos = MASS::ginv(B_rhos)

    C_i_rhos <- conf_normal_curvature(DELTA_pdispersao, Lobs, Lobs.aux = Lobs.aux_rhos)
    C_i_rhos_iden <- influential_index(C_i_rhos, C = 4)


    list_output[[which(perturbation=="dispersion")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                   C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                   C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }
}

if (any(perturbation == "response")){

  #calculo de las Deltas
  Delta_i <- lapply(N_i,
                    FUN = function(x){t(t(A%*%x) %*%
                                          (model$Dv-4*model$Dd) %*% A/phi)})
  Delta_phi<- t(A %*% model$residuals_mu) %*% ( model$Dv - 2* model$Dd) %*% A /(phi^2)
  Delta_rho <- lapply(B,
                      FUN = function(x){t( t(x%*%model$residuals_mu)%*%
                                             (4*model$Dd-model$Dv)%*% A -
                                             t(A%*%model$residuals_mu) %*%
                                             model$Dv %*% x) /phi})
  DELTA_pvarresposta <- t(cbind(do.call(cbind,Delta_i),
                                t(Delta_phi),
                                do.call(cbind,Delta_rho)))

  C_i <- conf_normal_curvature(DELTA_pvarresposta, Lobs, Lobs.aux = 0)
  C_i_iden <- influential_index(C_i, C = 4)

  list_output[[which(perturbation=="response")]] <- list()
  list_output[[which(perturbation=="response")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    # gamma

    B_gamma <- Lobs
    B_gamma[gamma_index,] <- 0
    B_gamma[,gamma_index] <- 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)

    C_i_gamma <- conf_normal_curvature(DELTA_pvarresposta, Lobs, Lobs.aux = Lobs.aux_gamma)
    C_i_gamma_iden <- influential_index(C_i_gamma, C = 4)

    # phi

    B_phi <- Lobs
    B_phi[phi_index,]<- 0
    B_phi[,phi_index]<- 0
    Lobs.aux_phi <- MASS::ginv(B_phi)

    C_i_phi <- conf_normal_curvature(DELTA_pvarresposta, Lobs, Lobs.aux = Lobs.aux_phi)
    C_i_phi_iden <- influential_index(C_i_phi, C = 4)

    # rhos

    B_rhos = Lobs
    B_rhos[rho_index,]<- 0
    B_rhos[,rho_index]<- 0
    Lobs.aux_rhos <- MASS::ginv(B_rhos)

    C_i_rhos <- conf_normal_curvature(DELTA_pvarresposta, Lobs, Lobs.aux = Lobs.aux_rhos)
    C_i_rhos_iden <- influential_index(C_i_rhos, C = 4)

    list_output[[which(perturbation=="response")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                  C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                  C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }

}

if (any(perturbation == "explanatory")){

  if (coef_index[1]==1){
    print("There is no explanatory variables")
    list_output[[which(perturbation=="explanatory")]] <- NULL
  } else {

    exp_index <- 2:coef_index[1]
    exp_numbers <- coef_index[1]-1

    list_output[[which(perturbation=="explanatory")]] <- list()

    for (r in exp_index){

      gamma_r<-model$f[[1]][r]

      #calculo de las Deltas
      Delta_i <- lapply(N_i,
                        FUN = function(x){t(t(A%*%x) %*%
                                              (4*model$Dd-model$Dv) %*% A/phi)})
      Delta_phi<- t(A %*% model$residuals_mu) %*% (2* model$Dd-model$Dv) %*% A /(phi^2)
      Delta_rho <- lapply(B,
                          FUN = function(x){t( t(x%*%model$residuals_mu)%*%
                                                 (model$Dv-4*model$Dd)%*% A -
                                                 #####
                                               t(A%*%model$residuals_mu) %*%
                                                 model$Dv %*% x) /phi})
      DELTA_pvarexpl <- gamma_r  * t(cbind(do.call(cbind,Delta_i),
                                           t(Delta_phi),
                                           do.call(cbind,Delta_rho)))

      C_i <- conf_normal_curvature(DELTA_pvarexpl, Lobs, Lobs.aux = 0)
      C_i_iden <- influential_index(C_i, C = 4)

      list_output[[which(perturbation=="explanatory")]][[r-1]] <- list()
      list_output[[which(perturbation=="explanatory")]][[r-1]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

      if(part){

        # gamma

        B_gamma <- Lobs
        B_gamma[gamma_index,] <- 0
        B_gamma[,gamma_index] <- 0
        Lobs.aux_gamma = MASS::ginv(B_gamma)

        C_i_gamma <- conf_normal_curvature(DELTA_pvarexpl, Lobs, Lobs.aux = Lobs.aux_gamma)
        C_i_gamma_iden <- influential_index(C_i_gamma, C = 4)

        # phi

        B_phi <- Lobs
        B_phi[phi_index,] <- 0
        B_phi[,phi_index] <- 0
        Lobs.aux_phi = MASS::ginv(B_phi)

        C_i_phi <- conf_normal_curvature(DELTA_pvarexpl, Lobs, Lobs.aux = Lobs.aux_phi)
        C_i_phi_iden <- influential_index(C_i_phi, C = 4)

        # rhos

        B_rhos <- Lobs
        B_rhos[rho_index,]<- 0
        B_rhos[,rho_index] <- 0
        Lobs.aux_rhos = MASS::ginv(B_rhos)

        C_i_rhos <- conf_normal_curvature(DELTA_pvarexpl,Lobs, Lobs.aux = Lobs.aux_rhos)
        C_i_rhos_iden <- influential_index(C_i_rhos, C = 4)

        list_output[[which(perturbation=="explanatory")]][[r-1]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                       C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                       C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)

      }
    }
  }

  names(list_output[[which(perturbation=="explanatory")]]) <- row.names(model$f[[1]])[-1]

}

names(list_output) <- perturbation

return(list_output)

}
