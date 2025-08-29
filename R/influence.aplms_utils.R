conf_normal_curvature <- function(DELTA, Lobs, Lobs.aux = 0){
  CC <- t(DELTA)%*%(solve(-Lobs)-Lobs.aux)%*%DELTA
  CC2 <- CC %*% CC
  CC_l <- CC/sqrt(tr(CC2))
  C_i <- diag(CC_l)
  return(C_i)
}


influential_index <- function(C_i){
  C_i_iden <-which(C_i>mean(C_i)+C*sd(C_i))
  return(C_i_iden)
}


