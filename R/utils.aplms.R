setupNiKi <- function(npc, basis, data, Knot) {
  k <- length(npc)
  ZZ <- list()
  N_i <- list()
  K_i <- list()
  for (i in 1:k) {
    XX <- as.list(substitute(list(npc[i])))[-1]
    YY <- s(XX, bs = basis[i], m = c(2, 3), k = Knot[i])
    YY$term <- npc[i]
    ZZ[[i]] <- mgcv::smoothCon(YY, data = data, absorb.cons = T)
    N_i[[i]] <- ZZ[[i]][[1]]$X
    K_i[[i]] <- (ZZ[[i]][[1]]$S)[[1]]
  }
  return(c(N_i, K_i))
}

backfitting <- function(N_i, A, Dv, k, f_init) {
  S_i <- list()
  S_i[[1]] <- tcrossprod(
    solve(t(A %*% N_i[[1]]) %*% Dv %*% (A %*% N_i[[1]])),
    (A %*% N_i[[1]])
  ) %*% Dv
  for (l in 1:k) {
    S_i[[l + 1]] <- tcrossprod(
      solve(t(A %*% N_i[[l + 1]]) %*% Dv %*% (A %*% N_i[[l + 1]]) + phi * lam[l] * K_i[[l + 1]]),
      (A %*% N_i[[l + 1]])
    ) %*% Dv
  }
  f <- f_init
  f0 <- S_i[[1]] %*% (A %*% (y - Reduce(`+`, mapply("%*%", N_i[-1], f[-1], SIMPLIFY = FALSE))))
  f[[1]] <- f0
  for (l in 1:k) {
    f_i <- S_i[[l + 1]] %*% (A %*% (y - Reduce(`+`, mapply("%*%", N_i[-(l + 1)], f[-(l + 1)], SIMPLIFY = FALSE))))
    f[[l + 1]] <- f_i
  }
  return(f)
}

pgam <- function(N_i, A, lam, phi, y, f_init) {
  AN_i <- lapply(N_i, FUN = function(x) {
    A %*% x
  })
  N_A <- do.call(cbind, AN_i)
  lam_list <- as.list(lam)
  lam_list <- append(0, lam_list)
  M_gamma <- mapply("*", lam_list, K_i, SIMPLIFY = FALSE)
  f_vector <- solve(t(N_A) %*% Dv %*% N_A + phi * Matrix::bdiag(M_gamma)) %*% t(N_A) %*% Dv %*% (A %*% y)
  f <- list()
  param_f <- sapply(f_init, nrow)
  f_dimension <- cumsum(param_f)
  f[[1]] <- cbind(f_vector[1:f_dimension[1]])
  for (l in 2:length(f_dimension)) {
    f[[l]] <- cbind(f_vector[(f_dimension[l - 1] + 1):f_dimension[l]])
  }
  return(f)
}