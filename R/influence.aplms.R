#' local influence of the object `aplms()`.
#'
#' Takes a fitted `aplms` object and outputs some diagnostic information about the fitting procedure and results. Returns the conformal normal curvature of the fitted `aplms` model object. The `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR` perturbations are available.
#'
#' @param model an object with the result of fitting additive partial linear models with symmetric errors.
#' @param perturbation
#' @param part a logical value to indicate whether the influential analysis is performed for $\gamma$, $\phi$ and $\rho$.
#' @param C The cutoff criterion such that $C_i > \bar{C_i} + C*sd(C_i)}$ to detect influential observations.
#' @param plot a logical value used to plot.
#' @import psych MASS
#'
influence.aplms(model,
                perturbation = c("case-weight","dispersion","response"), #,"explanatory", "corAR"),
                part = TRUE,
                C = 4,
                #ident = NULL,
                plot = FALSE,
                # labels = NULL,
                # iden = F
                ){

list_output <- list()

p = nrow(model$summary_table_phirho)-1
nn = nrow(model$data)
ONE<-cbind(rep(1,nn))
A <- matrix_A(model$summary_table_phirho[2:(p+1),1],nn)
B <- BB(p=4,nn)

phi <- model$summary_table_phirho[1,1]
N_i <- model$N_i

if (any(perturbation) == "case-weight"){

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
    Lobs <- model$LL_obs

    num_ger<-t(DELTA)%*%solve(-Lobs)%*%DELTA
    aux_ger <- num_ger %*% num_ger  # o considerar aux_ger<-SMUT::eigenMapMatMult(num_ger,num_ger)
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_ger/den_ger
    C_i <- diag(BB_ger)
    C_i_iden <-which(C_i>mean(C_i)+C*sd(C_i))

    list_output[[which(perturbation=="case-weight")]] <- list()
    list_output[[which(perturbation=="case-weight")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    coef_index<-c(sapply(N_i,dim)[2,],1,p)

    # gamma

    B_gamma = Lobs
    B_gamma[1:(sum(sapply(N_i,dim)[2,])),] = 0
    B_gamma[,1:(sum(sapply(N_i,dim)[2,]))] = 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)
    num_gamma<-t(DELTA)%*%(solve(-Lobs) -Lobs.aux_gamma)%*%DELTA
    aux_ger <- num_gamma %*% num_gamma # aux_ger <- SMUT::eigenMapMatMult(num_gamma,num_gamma)
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_gamma/den_ger
    C_i_gamma <- diag(BB_ger)
    C_i_gamma_iden <-which(C_i_gamma>mean(C_i_gamma)+C*sd(C_i_gamma))

    # phi

    B_phi = Lobs
    B_phi[(sum(coef_index)-p),]=0
    B_phi[,1:(sum(coef_index)-p)]=0
    Lobs.aux_phi = MASS::ginv(B_phi)

    num_phi<-t(DELTA)%*%(solve(-Lobs) -Lobs.aux_phi)%*%DELTA
    aux_ger <- num_phi %*% num_phi  ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_phi/den_ger
    C_i_phi <- diag(BB_ger)
    C_i_phi_iden <-which(C_i_phi>mean(C_i_phi)+C*sd(C_i_phi))


    # rhos

    B_rhos = Lobs
    B_rhos[((sum(coef_index)-p+1):sum(coef_index)),]=0
    B_rhos[,((sum(coef_index)-p+1):sum(coef_index))]=0
    Lobs.aux_rhos = MASS::ginv(B_rhos)

    num_rhos<-t(DELTA)%*%(solve(-Lobs) -Lobs.aux_rhos)%*%DELTA
    aux_ger <- num_rhos %*% num_rhos ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_rhos/den_ger
    C_i_rhos <- diag(BB_ger)
    C_i_rhos_iden<-which(C_i_rhos>mean(C_i_rhos)+4*sd(C_i_rhos))

    list_output[[which(perturbation=="case-weight")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                 C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                 C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }
}


if (any(perturbation) == "dispersion"){

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
  num_ger<-t(DELTA_pdispersao)%*%solve(-Lobs)%*%DELTA_pdispersao
  aux_ger <- num_ger %*% num_ger
  den_ger<-sqrt(tr(aux_ger))
  BB_ger<-num_ger/den_ger
  C_i <- diag(BB_ger)
  C_i_iden <-which(C_i>mean(C_i)+C*sd(C_i))

  list_output[[which(perturbation=="dispersion")]] <- list()
  list_output[[which(perturbation=="dispersion")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    coef_index<-c(sapply(N_i,dim)[2,],1,p)

    # gamma

    B_gamma = Lobs
    B_gamma[1:(sum(sapply(N_i,dim)[2,])),] = 0
    B_gamma[,1:(sum(sapply(N_i,dim)[2,]))] = 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)
    num_gamma<-t(DELTA_pdispersao)%*%(solve(-Lobs) -Lobs.aux_gamma)%*%DELTA_pdispersao
    aux_ger <- num_gamma %*% num_gamma # aux_ger <- SMUT::eigenMapMatMult(num_gamma,num_gamma)
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_gamma/den_ger
    C_i_gamma <- diag(BB_ger)
    C_i_gamma_iden <-which(C_i_gamma>mean(C_i_gamma)+C*sd(C_i_gamma))

    # phi

    B_phi = Lobs
    B_phi[(sum(coef_index)-p),]=0
    B_phi[,1:(sum(coef_index)-p)]=0
    Lobs.aux_phi = MASS::ginv(B_phi)

    num_phi<-t(DELTA_pdispersao)%*%(solve(-Lobs) -Lobs.aux_phi)%*%DELTA_pdispersao
    aux_ger <- num_phi %*% num_phi  ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_phi/den_ger
    C_i_phi <- diag(BB_ger)
    C_i_phi_iden <-which(C_i_phi>mean(C_i_phi)+C*sd(C_i_phi))


    # rhos

    B_rhos = Lobs
    B_rhos[((sum(coef_index)-p+1):sum(coef_index)),]=0
    B_rhos[,((sum(coef_index)-p+1):sum(coef_index))]=0
    Lobs.aux_rhos = MASS::ginv(B_rhos)

    num_rhos<-t(DELTA_pdispersao)%*%(solve(-Lobs) -Lobs.aux_rhos)%*%DELTA_pdispersao
    aux_ger <- num_rhos %*% num_rhos ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_rhos/den_ger
    C_i_rhos <- diag(BB_ger)
    C_i_rhos_iden<-which(C_i_rhos>mean(C_i_rhos)+4*sd(C_i_rhos))

    list_output[[which(perturbation=="dispersion")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                   C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                   C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }
}


if (any(perturbation) == "response"){

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

  num_ger<-t(DELTA_pvarresposta)%*%solve(-Lobs)%*%DELTA_pvarresposta
  aux_ger <- num_ger %*% num_ger  # o considerar aux_ger<-SMUT::eigenMapMatMult(num_ger,num_ger)
  den_ger<-sqrt(tr(aux_ger))
  BB_ger<-num_ger/den_ger
  C_i <- diag(BB_ger)
  C_i_iden <-which(C_i>mean(C_i)+C*sd(C_i))

  list_output[[which(perturbation=="response")]] <- list()
  list_output[[which(perturbation=="response")]][[1]] <- list(C_i=C_i, C_i_iden=C_i_iden)

  if(part){

    coef_index<-c(sapply(N_i,dim)[2,],1,p)

    # gamma

    B_gamma = Lobs
    B_gamma[1:(sum(sapply(N_i,dim)[2,])),] = 0
    B_gamma[,1:(sum(sapply(N_i,dim)[2,]))] = 0
    Lobs.aux_gamma = MASS::ginv(B_gamma)
    num_gamma<-t(DELTA_pvarresposta)%*%(solve(-Lobs) -Lobs.aux_gamma)%*%DELTA_pvarresposta
    aux_ger <- num_gamma %*% num_gamma # aux_ger <- SMUT::eigenMapMatMult(num_gamma,num_gamma)
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_gamma/den_ger
    C_i_gamma <- diag(BB_ger)
    C_i_gamma_iden <-which(C_i_gamma>mean(C_i_gamma)+C*sd(C_i_gamma))

    # phi

    B_phi = Lobs
    B_phi[(sum(coef_index)-p),]=0
    B_phi[,1:(sum(coef_index)-p)]=0
    Lobs.aux_phi = MASS::ginv(B_phi)

    num_phi<-t(DELTA_pvarresposta)%*%(solve(-Lobs) -Lobs.aux_phi)%*%DELTA_pvarresposta
    aux_ger <- num_phi %*% num_phi  ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_phi/den_ger
    C_i_phi <- diag(BB_ger)
    C_i_phi_iden <-which(C_i_phi>mean(C_i_phi)+C*sd(C_i_phi))


    # rhos

    B_rhos = Lobs
    B_rhos[((sum(coef_index)-p+1):sum(coef_index)),]=0
    B_rhos[,((sum(coef_index)-p+1):sum(coef_index))]=0
    Lobs.aux_rhos = MASS::ginv(B_rhos)

    num_rhos<-t(DELTA_pvarresposta)%*%(solve(-Lobs) -Lobs.aux_rhos)%*%DELTA_pvarresposta
    aux_ger <- num_rhos %*% num_rhos ## SMUT::eigenMapMatMult
    den_ger<-sqrt(tr(aux_ger))
    BB_ger<-num_rhos/den_ger
    C_i_rhos <- diag(BB_ger)
    C_i_rhos_iden<-which(C_i_rhos>mean(C_i_rhos)+4*sd(C_i_rhos))

    list_output[[which(perturbation=="response")]][[2]] <- list(C_i_gamma=C_i_gamma, C_i_gamma_iden=C_i_gamma_iden,
                                                                  C_i_phi=C_i_phi, C_i_phi_iden=C_i_phi_iden,
                                                                  C_i_rhos=C_i_rhos, C_i_rhos_iden=C_i_rhos_iden)
  }

}

names(list_output) <- perturbation

return(list_output)

}
