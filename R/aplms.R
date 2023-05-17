
aplms <- function(formula,npc=c("tdate","epi.week"), basis=c("cr","cc"),Knot=c(60,12),data,family=Normal,p=1,
                  control=list(tol=0.001,algorithm=c("BFGS")),
                  lam = c(100,10)){

  if (missingArg(formula)) {
    stop("The formula argument is missing.")
  }
  if (missingArg(npc)) {
    stop("The model needs at least one non-parametric component.")
  }
  if (missingArg(data)) {
    stop("The data argument is missing.")
  }
  if (missingArg(npc)) {
    stop("The model needs at least one non-parametric component.")
  }

  if(!all(npc %in% names(data))){
    stop("The non-parametric variables must be in data.")
  }

  if(control["algorithm"]=="BFGS"){
    gradient=NULL
  }
  if (control["algorithm"]=="Fisher.score"){
    gradient=grr
  }


  k <- length(npc)
  if (missingArg(basis)){
    if (k==1) basis <- c("cr")
    #if (k==2) basis <- c("cr","cc")
    if (k>=2) basis <- c("cr",rep("cc",k-1))
  }
  data1 <- model.frame(formula, data = data)

  y <- model.response(data1)

  N0 <- model.matrix(formula, data = data1)
  q <- ncol(N0)
  nn <- nrow(N0)

  K0<-matrix(0,nrow=ncol(N0),ncol=ncol(N0))


  ZZ<-mgcv::smoothCon(s(tdate,bs="cr",k=Knot[1]),data=data,absorb.cons=T)
  ZZ2<-mgcv::smoothCon(s(epi.week,bs="cc",k=Knot[2]),data=data,absorb.cons=T)

  N1<-ZZ[[1]]$X;K1<-(ZZ[[1]]$S)[[1]]
  N2<-ZZ2[[1]]$X;K2<-(ZZ2[[1]]$S)[[1]]

  N_i <- list(N0,N1,N2)
  K_i <- list(K0,K1,K2)

  rdf<-nrow(data)-ncol(N0)-p-1

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Symmetric error setup
  xi_t = family$g4(1, df = family$df,
                   alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                   sigmap = family$sigmap, k = family$k)
  # ww = family$g1(1, df = family$df,
  #                alpha = family$alpha, mp = family$mp, epsi = family$epsi,
  #                sigmap = family$sigmap, k = family$k)


  #initial values
  f0_ini =  rbind(mean(y),cbind(rep(0,dim(N0)[2]-1)))
  f1_ini =  cbind(rep(0,dim(N1)[2]))
  f2_ini =  cbind(seq(0,0,length=dim(N2)[2]))

  f_init = list(f0_ini,f1_ini,f2_ini)

  #rho_ini = c(0.52,0.2,0.1)
  rho_ini = rep(0.1,p)
  phi_ini = sd(y)/xi_t


  # if(is.null(lam)){ #lambda penalizadora
  #   #lam = rep(3550,length(npc))
  #   lam=c(100,10)
  # }

  #Estimation
  f0_aux = f0_ini
  f1_aux = f1_ini
  f2_aux = f2_ini

  f_aux<- list(f0_aux,f1_aux,f2_aux)
  #####
  conv_betaf=array()
  conv_betaf[1] = 1
  conv_geral=array()
  conv_geral[1] = 1
  i<-j<-1
  while(conv_geral[j] > 0.001){
    print(paste("While",j))
    i=1
    conv_betaf=array()
    conv_betaf[1] = 1
    j=j+1
    A = matrix_A(rho_ini,nn)
    while (conv_betaf[i] > 0.001){
      #a = res(y,  f0_ini,f1_ini, rho_ini,phi_ini,N_i)
      a = res(y, f_init, rho_ini, phi_ini,N_i)

      posicao = as.vector(family$g1(a, df = family$df,
                          alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                          sigmap = family$sigmap, k = family$k))
      Dv = (diag(-2*posicao))
      ############################
      print(i)
      i = i+1
      if(i>10){break;}############################control

      S_i <-list()
      S_i[[1]] = tcrossprod(solve(t(A%*% N_i[[1]] )%*% Dv %*%(A %*% N_i[[1]] )),
                            (A%*%N_i[[1]]))%*%Dv
      for(l in 1:k){
        S_i[[l+1]] = tcrossprod(solve(t(A%*%N_i[[l+1]])%*%Dv%*%(A%*%N_i[[l+1]])+phi_ini*lam[l]*K_i[[l+1]]),
                                (A%*%N_i[[l+1]]))%*%Dv
      }

      f = f_init
      f0 <- S_i[[1]]%*%(A%*%(y - Reduce(`+`, mapply('%*%',N_i[-1], f[-1], SIMPLIFY = FALSE))))
      f[[1]]<-f0
      for(l in 1:k){
        f_i <- S_i[[l+1]] %*%(A%*%(y - Reduce(`+`, mapply('%*%',N_i[-(l+1)], f[-(l+1)], SIMPLIFY = FALSE))))
        f[[l+1]]<- f_i
      }

      error <- mapply('-',f, f_init, SIMPLIFY = FALSE)
      conv_betaf[i] <- max(unlist(sapply(error,abs)))

      f_init = f
    }

    ##
    par1=optim(par=c(rho_ini,phi_ini),
               fn=logLik3.test,f=f_init,y=y,N_i=N_i,family=family, #gr=gradient,
               method ="L-BFGS-B",
               lower=c(rep(-1,p),0.001),upper=c(rep(1,p),Inf),
               control=list(fnscale=-1),hessian=T)

    dif_phi = par1$par[p+1] - phi_ini
    phi_ini = par1$par[p+1]

    dif_rho = par1$par[1:p] - rho_ini
    rho_ini = par1$par[1:p]


    f_error <- mapply('-',f_aux, f, SIMPLIFY = FALSE)
    f_aux <- f

    conv_geral[j] = max(abs(dif_phi),abs(dif_rho) ,unlist(sapply(f_error,abs)))

    if(j==25){break;}############################control
  }

  rho = rho_ini; phi = phi_ini
  A = matrix_A(rho,nn)
  a = res(y, f, rho, phi,N_i)
  posicao = as.vector(family$g1(a, df = family$df,
                      alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                      sigmap = family$sigmap, k = family$k))
  Dv = (diag(-2*posicao))
  Dm = diag(-2*posicao*as.vector(a)^2)

  Dd = diag( posicao*as.vector(a))

# loglik evaluation

  Lp <- logLik_fim.test(y, f, rho, phi, N_i, family)

  AN <- lapply(N_i, FUN = function(x) A %*% x)
  N_bar_a <- rlist::list.cbind(AN)                             #N_A
  N_bar<-rlist::list.cbind(N_i)
  K_ast<-phi*as.matrix(Matrix::bdiag(mapply('*',c(0,lam), K_i, SIMPLIFY = FALSE)))

  #
  H_alpha<-N_bar_a%*%solve((t(N_bar_a)%*%Dv)%*%N_bar_a+K_ast)%*%(t(N_bar_a)%*%Dv)        # (4.18)

  q1 <- t(N_bar_a)%*%Dv%*%N_bar_a
  dec <- eigen(q1)
  q12 <- dec$vectors%*%(diag(dec$values^(1/2)))%*%t(dec$vectors)
  q12m <- dec$vectors%*%(diag(dec$values^(-1/2)))%*%t(dec$vectors)
  auto <- q12m%*%(K_ast)%*%q12m

  df_alpha <- sum(1/(1+eigen(auto)$value))+p+1  #Efective degree of freedom

  effect <- solve(t(N_bar_a)%*%Dv%*%N_bar_a+K_ast)%*%t(N_bar_a)%*%Dv%*%N_bar_a   #pag 57, L4

  #effective degree of freedom per function
  n_i <- sapply(N_i,ncol)
  npc_dimension <- cumsum(n_i)
  dfk<-sum(diag(effect)[1:npc_dimension[1]])
  for(l in 2:length(npc_dimension)){
    dfk_i<-sum(diag(effect)[(npc_dimension[l-1]+1):npc_dimension[l]])
    dfk<-append(dfk,dfk_i)
  }

  ######
  II<-diag(nn)
  yhat<-H_alpha%*%(A%*%y)-(A-II)%*%y

  AIC<- -2*Lp+2*(df_alpha)
  BIC<- -2*Lp+log(nn)*(df_alpha)
  AICC<-AIC+1+
    2*(sum(diag(H_alpha))+1)/(nn-2-sum(diag(H_alpha)))
  AICc_alpha<-log((sum((sqrt(Dv)%*%(cbind(y-yhat)))^2))/nn)+1+
    2*(sum(diag(H_alpha))+1)/(nn-2-sum(diag(H_alpha)))


  muhat1<-Reduce(`+`, mapply('%*%',N_i, f, SIMPLIFY = FALSE))
  error_hat <- y-muhat1

  GCV<-((1/nn)*sum( ((sqrt(Dv)%*%(A%*%(y-muhat1)))^2)) )/
    ((1-sum(diag(H_alpha))/nn)^2)



  ######################
  dg_t <- family$g2(args, df = family$df, r = family$r,
                  s = family$s, alpha = family$alpha, mp = family$mp,
                  epsi = family$epsi, sigmap = family$sigmap,
                  k = family$k)
  fg_t <- family$g3(args, df = family$df, r = family$r,
                  s = family$s, alpha = family$alpha, mp = family$mp,
                  epsi = family$epsi, sigmap = family$sigmap,
                  k = family$k)
  ###################################################

  Dd1<-4*(dg_t/phi)*diag(nn)

  # Variancias de phi y rho.

  fis_phi = (nn/(4*phi^2))*(4*fg_t-1)       #pag 52
  #fis_rho = (4*dg_t*xi_t/(1-rho[1]^2))*(nn-1:p)

  VAR_phi<-1/fis_phi
  #VAR_rho<-1/fis_rho

  VAR_rho<-diag(solve(-par1$hessian))[1:p]

  # Variancias de las funciones suaves.

  const <- 4*(dg_t/phi)
  const2 <- Map("*",K_i,c(0,lam))

  rep_AN <- rep(AN,k+1)
  seq_AN <- rep(AN,each=k+1)

  fis_block <- mapply( FUN=function(X,Y) {t(X)%*% Y *const}, X= rep_AN , Y= seq_AN)

  diag_index <-diag(matrix(1:((k+1)^2),ncol=k+1))

  fis_block[diag_index] <- Map("+",fis_block[diag_index],const2)
  names(fis_block)<-letters[seq( from = 1, to = (k+1)^2 )]

  fis_FF <- list()
  for(l in 0:k){
    fis_FF[[l+1]] <- do.call(rbind,fis_block[(l*(k+1)+1):((l+1)*(k+1))])
  }

  fis_FF <- do.call(cbind,fis_FF)
  VAR_F<-solve(fis_FF)

  WALD <- list()
  for(l in 2:length(npc_dimension)){
    WALD[[l-1]]<- t(f[[l]]) %*%  solve(VAR_F[(npc_dimension[l-1]+1):npc_dimension[l],(npc_dimension[l-1]+1):npc_dimension[l]]) %*% f[[l]]
  }

  WALD_vec <- unlist(WALD)

  WALD_p_value <- mapply(FUN=function(X,Y){pchisq(q=X,df=Y,lower.tail = FALSE)},
                         X=WALD_vec,Y=dfk[-1])
#

  ###LL_observada

  c_i = as.vector(family$g5(a, df = family$df,
                            alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                            sigmap = family$sigmap, k = family$k))

  Dc = diag(c_i)
  Dd = diag( as.vector(a) * c_i)

  LL_block <- mapply( FUN=function(X,Y) {(1/phi)*t(X)%*% (4*Dd-Dv) %*% Y }, X= rep_AN , Y= seq_AN)

  diag_index <-diag(matrix(1:((k+1)^2),ncol=k+1))

  LL_block[diag_index] <- Map("-",LL_block[diag_index],const2)
  names(LL_block)<-letters[seq( from = 1, to = (k+1)^2 )]

  LL_FF <- list()
  for(l in 0:k){
    LL_FF[[l+1]] <- do.call(rbind,LL_block[(l*(k+1)+1):((l+1)*(k+1))])
  }

  LL_FF <- do.call(cbind,LL_FF)


  delta_i <- a^2
  ONE <- cbind(rep(1,nn))

  LL_phi <- (1/phi^2)*(nn/2 + t(delta_i)%*%Dc%*%delta_i - t(delta_i)%*%Dv%*%ONE)

  B <- BB(p=p,nn=nn)
  rep_B <- rep(B,p)
  seq_B <- rep(B,each=p)
  LL_rho_block <- mapply( FUN=function(X,Y){
    (1/phi)*(t(X%*%error_hat) %*% (4*Dd-Dv) %*% (Y%*%error_hat))
  }, X= rep_B , Y= seq_B)

  LL_rho<- matrix(LL_rho_block,nrow=p)

  LL_FF_phi_block <- lapply(AN, FUN=function(x) {(1/phi^2)* (t(x)%*% (2*Dd-Dv) %*% (A%*% error_hat))})

  LL_FF_phi <- do.call(rbind,LL_FF_phi_block)

  rep_N_i <- rep(N_i,each=p)
  seq_B <- rep(B,k+1)
  LL_FF_rho_block <- mapply( FUN=function(X,Y) {
    (1/phi)*((t( A%*% X)%*% (Dv-4*Dd) %*% (Y %*%error_hat))+ t(Y %*% X) %*%  Dv %*% A %*% error_hat)
  },
  X= rep_N_i , Y= seq_B)

  LL_FF_rho <- list()
  for(l in 0:k){
    LL_FF_rho[[l+1]] <- do.call(cbind,LL_FF_rho_block[(l*(p)+1):((l+1)*p)])
  }

  LL_FF_rho <- do.call(rbind,LL_FF_rho)


  LL_phi_rho_block <- lapply(B, FUN=function(x) {(1/phi^2)* (t(x%*% error_hat)%*% (Dv-2*Dd) %*% (A%*% error_hat))})
  LL_phi_rho <- do.call(cbind,LL_phi_rho_block)

  LL_obs<- rbind(
    cbind(LL_FF,LL_FF_phi,LL_FF_rho),
    cbind(t(LL_FF_phi),LL_phi,LL_phi_rho),
    cbind(t(LL_FF_rho),t(LL_phi_rho),LL_rho))


  #Residuals

  #Pearson Residuals
  res_per <- (y - muhat1)/(sqrt(xi_t*phi))
  #qqnorm(res_per);qqline(res_per)

  #Quantile Residuals           #ver ver esto para distribuciones diferentes.
  #res_quant<-c(qnorm(rmutil::ppowexp((y-yhat),m=0,s=phi1,f=1/(1+kk))))

  #coeficientes

  # concise<-FALSE
  # digits = max(3L, getOption("digits") - 3L)

  est_coef <- as.vector(f0)
  ee <- sqrt(diag(VAR_F)[1:length(f0)])
  t_test <- est_coef/ee
  # p_value <- format.pval(sapply(t_test, FUN = function(x){2 * pt(abs(x), rdf, lower.tail = FALSE)}),
  #                        digits = digits, if (!concise) .Machine$double.eps else 1e-4)
  # summary_table<-cbind(sprintf('%.6f',est_coef),
  #                    sprintf('%.6f',ee),
  #                    sprintf('%.6f',t_test),
  #                    p_value)

  p_value <- sapply(t_test, FUN = function(x){2 * pt(abs(x), rdf, lower.tail = FALSE)})

  summary_table<-cbind(est_coef,
                       ee,
                       t_test,
                       p_value)

  terms_formula<-stats::terms(formula)
  var_names <- attr(terms_formula, "term.labels")


  rownames(summary_table)<- c("intercept",var_names)
  colnames(summary_table) <- c("Estimate","Std. Error","t value","Pr(>|t|)")

  WALD_f <- cbind(dfk[-1],WALD_vec, WALD_p_value)
  rownames(WALD_f)<- npc
  colnames(WALD_f)<- c("df","Wald","p-value")

  WALD_rho <- rho/sqrt(VAR_rho)

  # p_value_rho_normal <- format.pval(sapply(WALD_rho, FUN = function(x){2 * pnorm(abs(x), lower.tail = FALSE)}),
  #                        digits = digits, if (!concise) .Machine$double.eps else 1e-4)
  # p_value_rho_t <- format.pval(sapply(WALD_rho, FUN = function(x){2 * pt(x, rdf, lower.tail = FALSE)}),
  #                                   digits = digits, if (!concise) .Machine$double.eps else 1e-4)
  # summary_table_rho <- cbind(sprintf('%.6f',rho),
  #                            sprintf('%.6f',sqrt(VAR_rho)),
  #                            sprintf('%.6f',WALD_rho),
  #                            p_value_rho_normal,
  #                            p_value_rho_t)

  p_value_rho_normal <- sapply(WALD_rho, FUN = function(x){2 * pnorm(abs(x), lower.tail = FALSE)})
  p_value_rho_t <- sapply(WALD_rho, FUN = function(x){2 * pt(x, rdf, lower.tail = FALSE)})
  summary_table_rho <- cbind(rho,
                             sqrt(VAR_rho),
                             WALD_rho,
                             p_value_rho_normal,
                             p_value_rho_t)


  colnames(summary_table_rho)<- c("rho","ee","Wald","p-value_normal","p-value_t")


  WALD_phi <- phi/sqrt(VAR_phi)

  # p_value_phi_normal <- format.pval(2*pnorm(abs(WALD_phi),lower.tail = FALSE),
  #                                   digits = digits, if (!concise) .Machine$double.eps else 1e-4)
  # p_value_phi_t <- format.pval(2*pt(abs(WALD_phi),rdf,lower.tail = FALSE),
  #                              digits = digits, if (!concise) .Machine$double.eps else 1e-4)
  # summary_table_phi <- c(sprintf('%.6f',phi),
  #                            sprintf('%.6f',sqrt(VAR_phi)),
  #                            sprintf('%.6f',WALD_phi),
  #                            p_value_phi_normal,
  #                            p_value_phi_t)

  p_value_phi_normal <- 2*pnorm(abs(WALD_phi),lower.tail = FALSE)
  p_value_phi_t <- 2*pt(abs(WALD_phi),rdf,lower.tail = FALSE)
  summary_table_phi <- c(phi,
                         sqrt(VAR_phi),
                         WALD_phi,
                         p_value_phi_normal,
                         p_value_phi_t)

  summary_table_rhophi <- rbind(summary_table_rho,summary_table_phi)

  rownames(summary_table_rhophi)<- c(paste0("rho",as.character(1:p)),"phi")

  fit<-list(
    formula = formula,
    family = family,
    npc = npc,
    lam = lam,
    coef= summary_table,
    VAR_F = VAR_F,
    WALD_f = WALD_f,
    parametric_error =summary_table_rhophi,
    N_i = N_i,
    f= f,
    Dv = Dv,
    Dm = Dm,
    Dc = Dc,
    Dd = Dd,
    delta = delta_i,
    LL_obs = LL_obs,
    loglike = Lp,
    total_df = df_alpha,
    parametric_df = dfk[1],
    npc_df = dfk[-1],
    AIC = AIC,
    BIC = BIC,
    AICC = AICC,
    AICc_alpha = AICc_alpha,
    GCV = GCV,
    yhat = yhat,
    muhat = muhat1,
    residuals_y = y-yhat,
    residuals_mu = y-muhat1
  )
  class(fit) <- "aplms"
  invisible(fit)
}



print.aplms<- function(x,...){
  if(class(x) != "aplms")
    stop("not a aplms object")
  cat("Call:\n")
  print(x$this.call)

}


