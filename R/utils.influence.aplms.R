#Delta computation depending on the perturbation scheme

influence_DELTA <- function(model,
                            perturb_scheme = "case-weight",
                            r = 1,  # index of explanatory variable
                            k = 1   # index of autocorrelation coefficient
                              ){
  p <- nrow(model$summary_table_phirho)-1
  nn <- nrow(model$data)
  ONE <- cbind(rep(1,nn))
  A <- matrix_A(model$summary_table_phirho[2:(p+1),1],nn)
  B <- BB(p=p,nn)

  phi <- model$summary_table_phirho[1,1]
  N_i <- model$N_i
  Lobs <- model$LL_obs

  if(perturb_scheme == "case-weight"){
    Delta_i <- lapply(model$N_i,
                      FUN = function(x){
                        t(t(A%*%x) %*% model$Dv %*% diag(c(A %*% model$residuals_mu))/phi)
                      })
    Delta_phi<- -ONE/(2*phi) + model$Dm%*%ONE/(2*phi)
    Delta_rho <- lapply(B,
                        FUN = function(x){t(t(x %*% model$residuals_mu) %*%
                                              model$Dv %*% diag(c(A %*% model$residuals_mu))/phi)})
  }


  if(perturb_scheme == "dispersion"){
    Delta_i <- lapply(N_i,
                      FUN = function(x){t(t(A%*%x) %*% (model$Dv-2*model$Dd) %*%
                                            diag(c(A %*% model$residuals_mu))/phi)})
    Delta_phi<- t(1/(2*phi)*t(model$delta) %*% ( model$Dv - 2* model$Dd))
    Delta_rho <- lapply(B,
                        FUN = function(x){t( t(x%*%model$residuals_mu)%*%
                                               (2*model$Dd-model$Dv)%*%
                                               diag(c(A%*%model$residuals_mu)) /phi)})
  }


  if(perturb_scheme == "response"){
    Delta_i <- lapply(N_i,
                      FUN = function(x){t(t(A%*%x) %*%
                                            (model$Dv-4*model$Dd) %*% A/phi)})
    Delta_phi<- t(t(A %*% model$residuals_mu) %*% ( model$Dv - 2* model$Dd) %*% A /(phi^2))
    Delta_rho <- lapply(B,
                        FUN = function(x){t( t(x%*%model$residuals_mu)%*%
                                               (4*model$Dd-model$Dv)%*% A -
                                               t(A%*%model$residuals_mu) %*%
                                               model$Dv %*% x) /phi})
  }

  if(perturb_scheme == "explanatory"){

    gamma_r<-model$f[[1]][r+1]

    Delta_i <- lapply(N_i,
                      FUN = function(x){gamma_r *t(t(A%*%x) %*%
                                            (4*model$Dd-model$Dv) %*% A/phi)})
    Delta_phi<- t(gamma_r  * t(A %*% model$residuals_mu) %*% (2* model$Dd-model$Dv) %*% A /(phi^2))
    Delta_rho <- lapply(B,
                        FUN = function(x){gamma_r  * t( t(x%*%model$residuals_mu)%*%
                                               (model$Dv-4*model$Dd)%*% A -
                                               #####
                                             t(A%*%model$residuals_mu) %*%
                                               model$Dv %*% x) /phi})
  }


  if(perturb_scheme == "corAR"){
    Br <- B[[k]]
    Delta_i <- lapply(N_i,
                      FUN = function(x){(1/phi)*t(t(A%*%x) %*%
                                                    (model$Dv-4*model$Dd) %*% diag(c(Br%*%model$residuals_mu)) +
                                                    t(Br%*%x) %*% model$Dv %*% diag(c(A%*%model$residuals_mu)))})

    Delta_phi<- t(t(A %*% model$residuals_mu) %*% ( model$Dv - 2* model$Dd) %*% (diag(c(Br%*%model$residuals_mu))) /(phi^2))

    Delta_rho <- lapply(B,
                        FUN = function(x){(1/phi) *t( t(x%*%model$residuals_mu)%*%
                                                        (4*model$Dd-model$Dv)%*%
                                                        (diag(c(Br%*%model$residuals_mu))))})
  }

  DELTA <- t(cbind(do.call(cbind,Delta_i),
                   Delta_phi,
                   do.call(cbind,Delta_rho)))

  return(DELTA)
}


set_zero_matrix <- function(mat, rows = NULL, cols = NULL) {

  if (!is.matrix(mat)) stop("Input must be a matrix.")
  if (!is.null(rows) && any(rows > nrow(mat) | rows < 1)) stop("Row indices out of range.")
  if (!is.null(cols) && any(cols > ncol(mat) | cols < 1)) stop("Column indices out of range.")
  if (!is.null(rows)) {
    mat[rows, ] <- 0
  }
  if (!is.null(cols)) {
    mat[, cols] <- 0
  }
  return(mat)
}


Gp <- function(Lobs, index = NULL){

  Lobs <- set_zero_matrix(mat = Lobs,
                             rows = index,
                             cols = index)
  Lobs.aux = MASS::ginv(Lobs)
  return(Lobs.aux)
}





conf_normal_curvature <- function(DELTA, Lobs, Lobs.aux = 0){
  CC <- t(DELTA)%*%(solve(-Lobs)-Lobs.aux)%*%DELTA
  CC2 <- CC %*% CC
  CC_l <- CC/sqrt(sum(diag(CC2)))
  C_i <- diag(CC_l)
  return(C_i)
}


influential_index <- function(C_i, C = 4){
  C_i_iden <-which(C_i > mean(C_i, na.rm=  TRUE) + C * sd(C_i, na.rm = TRUE))
  return(C_i_iden)
}


