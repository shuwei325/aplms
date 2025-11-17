# Utility Functions to support the function `influence.aplms`

# Delta computation depending on the perturbation scheme
#
# @param model an object with the result of fitting additive partial linear models with symmetric errors.
# @param perturb_scheme A string vector specifying a perturbation scheme: `case-weight`, `dispersion`, `response`, `explanatory`, and `corAR`.
# @param r Index of explanatory variable
# @param k Index of autocorrelation coefficient
# @return A matrix Delta to support the perturbation scheme calculation.
influence_DELTA <- function(model,
                            perturb_scheme = "case-weight",
                            r = 1,
                            k = 1
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

# Setting the columns and rows indices of a matrix to zero
# @param mat squared matrix
# @param rows Row index
# @param cols Column Index
# @return A matrix with the specified rows and columns set to zero.
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

# Setting columns and rows indices of a matrix to zero and compute the Generalized Inversion.
# @param Lobs The observed information matrix
# @param Column and row indices to be set to zero.
# @return The generalized inverse of the matrix \code{Lobs}, with the specified rows and columns set to zero, to support perturbation scheme computations.
Gp <- function(Lobs, index = NULL){

  Lobs <- set_zero_matrix(mat = Lobs,
                             rows = index,
                             cols = index)
  Lobs.aux = MASS::ginv(Lobs)
  return(Lobs.aux)
}


# Conformal normal curvature of local influence
#
# Computes the conformal normal curvature.
#
# @param DELTA A numeric matrix, \eqn{\Delta}, that depends on the perturbation scheme.
# @param Lobs A square numeric matrix representing the observed information matrix.
# @param Lobs.aux An optional numeric matrix of the same dimension as
#   \code{Lobs}. Default is \code{0} (no auxiliary adjustment).
# @return A vector, results of conformal normal curvature.
conf_normal_curvature <- function(DELTA, Lobs, Lobs.aux = 0){
  CC <- t(DELTA)%*%(solve(-Lobs)-Lobs.aux)%*%DELTA
  CC2 <- CC %*% CC
  CC_l <- CC/sqrt(sum(diag(CC2)))
  C_i <- diag(CC_l)
  return(C_i)
}


influential_plot1 <- function(k, output, labels = NULL, C = 4,...){
  C_i <- output[[k]]
  name <- names(output)[[k]]
  C_mean <- mean(C_i, na.rm = TRUE)
  C_sd <- sd(C_i, na.rm = TRUE)
  thres <- C_mean + C * C_sd
  df <- data.frame(labels = labels, C_i = C_i)
  df$label_text <- ifelse(C_i > thres, labels, NA)
  p <- ggplot(df, aes(x = labels, y = C_i)) +
    geom_point() +
    geom_hline(yintercept = thres, color = "red") +
    geom_text(aes(label = label_text), vjust = -0.5, na.rm = TRUE) +
    coord_cartesian(ylim = range(0, max(C_i)+C_sd)) +
    labs(x = "t", y = name) +
    theme_minimal()
  return(p)
}
