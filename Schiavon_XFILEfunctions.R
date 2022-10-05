#---------------------------------------------------------------------------------------------#
# ---- XFILE ALGORITHM R FUNCTIONS ----
# Copyright: Lorenzo Schiavon 2021,
# Licensed under Creative Commons Attribution-NonCommercial 4.0
#---------------------------------------------------------------------------------------------#


sigmaNullFunctions=function(a_sigma, b_sigma){
  
  WUpdate = function(reps, predictor1, param, predictor2, eps){
    out = (pmax(predictor1,0)+eps)*rep((pmax(predictor2,0)+eps)*param, each=reps)
    return(out)
  }
  D2invUpdate = function(reps, W, barz, param0, b_sigma){
    out = W^2/(2*b_sigma + (barz -rep(param0,reps))^2)
    return( out )
  }
  WUpdate_beta = function(reps, predictor, param1, param2, eps){
    out = param1*rep((pmax(predictor,0)+eps)*param2, each=reps)
    return(out)
  }
  D2invUpdate_beta = function(reps, W, barz, predictor0, eps){
    out = W^2/(2*b_sigma + (barz - rep( (pmax(predictor0,0)+eps),reps))^2)
    return(out)
  }
  constant = (a_sigma +0.5)
  variance_constant = b_sigma
  
  loglikelihood = function(data, mu){
    -constant*log(1+(data-mu)^2/(2*variance_constant))
  }
  
  return(list(WUpdate=WUpdate, D2invUpdate = D2invUpdate,
              WUpdate_beta=WUpdate_beta, D2invUpdate_beta=D2invUpdate_beta, constant=constant,
              loglikelihood=loglikelihood))

}


# every step reduce the logposterior but for the Beta step
# we accept the new beta only if they reduce the logposterior
xfileCoordinateMM = function(y, n=nrow(y), p=ncol(y), x=NULL, w=NULL, M_h=matrix(0,n,p), hyperparameter, 
                             data_distribution ="gaussian", 
                             listfunction = sigmaNullFunctions(hyperparameter$a_sigma, hyperparameter$b_sigma),
                             tolerance = 10^(-5), max_it = 100, seed=1111){
  
  require("invgamma")
  if(is.null(x)){x = matrix(1,n,1)}
  if(is.null(w)){w = matrix(1,p,1)}
  qx = ncol(x)
  qw = ncol(w)

  if(data_distribution %in% c("gaussian", "Gaussian", "normal", "Normal")){
    y=y-M_h
    wy0 = NULL
    zupdate = function(y, M_h, F_h, wy0){return(y)}
  }else if(data_distribution %in% c("binary", "Binary", "probit", "Probit")){
    wy0 = which(y==0)
    zupdate = function(y, M_h, F_h, wy0){
      z= M_h
      all = length(z)
      FM = which(F_h< -M_h)
      wF_h = which((FM %in% wy0)|( all[-FM] %in% all[-wy0] ))
      if(length(wF_h)>0){
        z[wF_h] = F_h[wF_h]  
      }
      return(z)
    }
    Sigma = diag(n*p)
  }else if(data_distribution %in% c("censored", "Censored", "Truncated", "Truncated")){
    wy0 = which(y==0)
    y[-wy0] = y[-wy0]-M_h[-wy0]
    zupdate = function(y, M_h, F_h, wy0){
      z=y
      z[wy0] = -M_h[wy0]
      wF_h = wy0[which(F_h[wy0]< -M_h[wy0])]
      if(length(wF_h)>0){
        z[wF_h] = F_h[wF_h]  
      }
      return(z)
    }
    
  }else{stop("The argument of data_distribution must be either 'gaussian', 'binary' or 'censored'")}

  # useful stuff
  lratio_cn = log(hyperparameter$c_n)-log(1-hyperparameter$c_n)
  lratio_cp = log(hyperparameter$c_p)-log(1-hyperparameter$c_p)
  eps = hyperparameter$epsilon
  meanBetaeta = rep(0,qx); meanBetaeta[1]= 1-eps
  meanBetalambda = rep(0,qw); meanBetalambda[1]= 1-eps
  constant_update = 0.5/listfunction$constant
  
  # save tracking
  Vartheta = rep(NA, max_it+1)
  Z = array(NA, dim = c(max_it+1, n, p))
  Eta = PhiEta = matrix(NA, max_it+1, n)
  Lambda = PhiLambda = matrix(NA, max_it+1, p)
  BetaEta = matrix(NA, max_it+1, qx)
  BetaLambda = matrix(NA, max_it+1,qw)
  lpost = rep(NA, max_it+1)
  
  # first z
  F_h = matrix(0,n,p)
  Z[1,,] = z = zupdate(y, M_h, F_h, wy0)
  loglikz = listfunction$loglikelihood(z, mu=0)
  
  # we want to converge to the posterior mode with highest prior
  # initialization close to prior mode 
  set.seed(seed)
  Vartheta[1] = vartheta = hyperparameter$b_theta/(hyperparameter$a_theta+1)
  Eta[1,] = eta = rnorm(n,0, 1)
  PhiEta[1,] = phieta = rep(1,n)
  BetaEta[1,] = betaeta = meanBetaeta
  lambda = rnorm(p,0,sqrt(vartheta))
  Lambda[1,] = lambda/sqrt(vartheta)
  PhiLambda[1,] = philambda = rep(1, p)
  BetaLambda[1,] = betalambda = meanBetalambda
  predeta = drop(x%*%betaeta)
  predlambda = drop(w%*%betalambda)
  wphieta = seq_len(n)
  wphilambda = seq_len(p)
  lwphieta = n
  lwphilambda = p
  
  
  # first lpost
  F_h = matrix((pmax(predeta,0)+eps)*eta*phieta*
           rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
  lpost[1] = sum(listfunction$loglikelihood(z, mu = F_h), na.rm = T) +
    dinvgamma(vartheta, hyperparameter$a_theta, hyperparameter$b_theta,log=T)+
    sum(dnorm(eta,0,1, log=T)) + sum(dnorm(lambda,0,sqrt(vartheta), log=T)) +
    sum(dnorm(betaeta, meanBetaeta, 1, log=T)) + sum(dnorm(betalambda, meanBetalambda,1, log=T)) +
    sum(log(phieta*hyperparameter$c_n + (1-phieta)*hyperparameter$c_n))+
    sum(log(philambda*hyperparameter$c_p + (1-philambda)*hyperparameter$c_p))
  
  t=1; lp_decrease = tolerance+1
  while((t<max_it+1)&(lp_decrease>tolerance)){
    
    # update eta
    if(lwphilambda>0){
      Weta = listfunction$WUpdate(n, predeta, lambda[wphilambda], predlambda[wphilambda], eps)
      barz = c(z[,wphilambda])/Weta
      D2inveta = listfunction$D2invUpdate(lwphilambda, Weta, barz, eta, hyperparameter$b_sigma)
      Diageta = (rowSums(matrix(D2inveta,n,lwphilambda), na.rm = T)+ constant_update)^{-1}
      eta = Diageta*rowSums(matrix(D2inveta*barz, n,lwphilambda), na.rm=T)
    }else{ eta = rnorm(n,0, 0.1)}
          
    # update phieta
    Fh_phieta1 = matrix((pmax(predeta,0)+eps)*eta*
                          rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
    lphieta = rowSums(loglikz-listfunction$loglikelihood(z, mu=Fh_phieta1), na.rm =T)
    wphieta.tmp = which(lphieta < lratio_cn)
    lwphieta.tmp = length(wphieta.tmp)
    # if every row is equal to zero: we skip the iteration but removing lambda = 0
    if(lwphieta.tmp>0){
      wphieta = wphieta.tmp; lwphieta = lwphieta.tmp
    }else{
      wphieta = setdiff(wphieta, which(eta==0))
      lwphieta = length(wphieta)
    }
    phieta = rep(0,n)
    phieta[wphieta]=1
  
    # update betaeta
    Seta = which((predeta>0)&(phieta>0))
    if((lwphilambda>0)&(length(Seta)>0)){
      Wbetaeta = listfunction$WUpdate_beta(length(Seta), predlambda[wphilambda], eta[Seta], lambda[wphilambda], eps)
      barzseta = c(z[Seta,wphilambda])/Wbetaeta
      D2invbetaeta = listfunction$D2invUpdate_beta(lwphilambda, Wbetaeta, barzseta, predeta[Seta], eps)
      xD2inv = D2invbetaeta *matrix(rep(t(x[Seta,]), lwphilambda), ncol=qx, byrow=TRUE)
      xD2inv[is.na(xD2inv)] = 0  # to deal with the NA in crossprod
      Diagetabeta = solve(crossprod(matrix(rep(t(x[Seta,]), lwphilambda), ncol=qx,byrow=TRUE), xD2inv)+ diag(constant_update,qx))
      barzseta[is.na(barzseta)] = 0  # to deal with the NA in crossprod
      betaeta.tmp = Diagetabeta%*%(crossprod(xD2inv, barzseta)+ constant_update*meanBetaeta)
      if(!is.nan(sum(betaeta.tmp))){
        predeta.tmp = drop(x%*%betaeta.tmp)
        F_h.tmp = matrix((pmax(predeta.tmp,0)+eps)*eta*phieta*
                           rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
        F_h = matrix((pmax(predeta,0)+eps)*eta*phieta*
                       rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
        if(sum(listfunction$loglikelihood(z, mu = F_h), na.rm = T)+ sum(dnorm(betaeta, meanBetaeta,1, log=T))<
           sum(listfunction$loglikelihood(z, mu = F_h.tmp), na.rm = T)+ sum(dnorm(betaeta.tmp, meanBetaeta,1, log=T))){
          betaeta = betaeta.tmp
        }else{
          betaeta = betaeta + (betaeta.tmp-betaeta)/10
        }
      }
    }else{
      betaeta = meanBetaeta
    }
    predeta = drop(x%*%betaeta)
    
    # update lambda
    if(lwphieta>0){
      Wlambda = listfunction$WUpdate(p, predlambda, eta[wphieta], predeta[wphieta], eps)
      barzt = c(t(z[wphieta,]))/Wlambda
      D2invlambda = listfunction$D2invUpdate(lwphieta, Wlambda, barzt, lambda, hyperparameter$b_sigma)
      Diaglambda = (colSums(matrix(D2invlambda,lwphieta, p, byrow=T), na.rm = T) + constant_update/vartheta)^{-1}
      lambda = Diaglambda*colSums(matrix(D2invlambda*barzt, lwphieta ,p, byrow = T), na.rm = T)
    }else{
      lambda = rnorm(p, 0,0.01)
    }
    
    # update philambda. 
    Fh_philambda1 = matrix((pmax(predeta,0)+eps)*eta*phieta*
                          rep((pmax(predlambda,0)+eps)*lambda, each=n), n, p)
    lphilambda = colSums(loglikz-listfunction$loglikelihood(z, mu=Fh_philambda1), na.rm = T)
    wphilambda.tmp = which(lphilambda<lratio_cp )
    lwphilambda.tmp = length(wphilambda.tmp)
    # if every col is equal to zero: we skip the iteration but removing lambda = 0
    if(lwphilambda.tmp>0){
      wphilambda = wphilambda.tmp; lwphilambda = lwphilambda.tmp
    }else{
      wphilambda = setdiff(wphilambda, which(lambda==0))
      lwphilambda = length(wphilambda)
    }
    philambda = rep(0,p)
    philambda[wphilambda]=1 
    
    # update betalambda
    Slambda = which((predlambda>0)&(philambda>0))
    if((lwphieta>0)&(length(Slambda)>0)){
      Wbetalambda = listfunction$WUpdate_beta(length(Slambda), predeta[wphieta], lambda[Slambda], eta[wphieta], eps)
      tbarzslambda = c(t(z[wphieta, Slambda]))/Wbetalambda
      D2invbetalambda = listfunction$D2invUpdate_beta(lwphieta, Wbetalambda, tbarzslambda, predlambda[Slambda], eps)
      wD2inv = D2invbetalambda*matrix(rep(t(w[Slambda,]), lwphieta),ncol=qw,byrow=TRUE)
      wD2inv[is.na(wD2inv)] = 0  # to deal with the NA in crossprod
      Diaglambdabeta = solve(crossprod(matrix(rep(t(w[Slambda,]), lwphieta),ncol=qw,byrow=TRUE),wD2inv)+ diag(constant_update,qw))
      tbarzslambda[is.na(tbarzslambda)] = 0  # to deal with the NA in crossprod
      betalambda.tmp = Diaglambdabeta%*%(crossprod(wD2inv, tbarzslambda)+ constant_update*meanBetalambda)
      if(!is.nan(sum(betalambda.tmp))){
        predlambda.tmp = drop(w%*%betalambda.tmp)
        F_h.tmp = matrix((pmax(predeta,0)+eps)*eta*phieta*
                           rep((pmax(predlambda.tmp,0)+eps)*lambda*philambda, each=n), n, p)
        F_h = matrix((pmax(predeta,0)+eps)*eta*phieta*
                       rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
        if(sum(listfunction$loglikelihood(z, mu = F_h), na.rm = T)+ sum(dnorm(betalambda, meanBetalambda,1, log=T))<
           sum(listfunction$loglikelihood(z, mu = F_h.tmp), na.rm = T)+ sum(dnorm(betalambda.tmp, meanBetalambda,1, log=T))){
          betalambda = betalambda.tmp
        }else{
          betalambda = betalambda + (betalambda.tmp-betalambda)/10
        }
      }
    }else{
      betalambda = meanBetalambda
    }
    predlambda = drop(w%*%betalambda)
    
    # update vartheta
    vartheta = (hyperparameter$b_theta + 0.5*sum(lambda^2))/(hyperparameter$a_theta+0.5*p+1)
    
    # update factor mean
    F_h = matrix((pmax(predeta,0)+eps)*eta*phieta*
      rep((pmax(predlambda,0)+eps)*lambda*philambda, each=n), n, p)
    
    # update z
    z = zupdate(y, M_h, F_h, wy0 = wy0)
    loglikz = listfunction$loglikelihood(z, mu=0)
    
    # save
    Vartheta[t+1] = vartheta
    Z[t+1,,] = z
    Eta[t+1,] = eta
    PhiEta[t+1,] = phieta
    BetaEta[t+1,] = betaeta
    Lambda[t+1,] = lambda/sqrt(vartheta)
    PhiLambda[t+1,] = philambda
    BetaLambda[t+1,] = betalambda
    
    # log posterior
    lpost[t+1] = sum(listfunction$loglikelihood(z, mu = F_h), na.rm =T) +
      dinvgamma(vartheta, hyperparameter$a_theta, hyperparameter$b_theta,log=T)+
      sum(dnorm(eta,0,1, log=T)) + sum(dnorm(lambda,0,sqrt(vartheta), log=T)) +
      sum(dnorm(betaeta, meanBetaeta, 1, log=T)) + sum(dnorm(betalambda, meanBetalambda,1, log=T)) +
      sum(log(phieta*hyperparameter$c_n + (1-phieta)*hyperparameter$c_n))+
      sum(log(philambda*hyperparameter$c_p + (1-philambda)*hyperparameter$c_p))
    
    # log posterior decrease ratio
    lp_decrease = abs((lpost[t+1]-lpost[t])/lpost[t])

    t=t+1
  }
  
  if(t<max_it+1){
    Vartheta = Vartheta[1:t]
    Z = Z[1:t,,]
    Eta = Eta[1:t,]
    PhiEta = PhiEta[1:t,]
    BetaEta[1:t,]
    Lambda[1:t,]
    PhiLambda[1:t,]
    BetaLambda[1:t,]
    lpost = lpost[1:t]
  }else{warning("Maximum iteration reached: algorithm stops")}
  
  return(list(F_h = F_h, Vartheta = Vartheta, Z = Z, Eta = Eta, PhiEta = PhiEta, BetaEta = BetaEta,
              Lambda = Lambda, PhiLambda = PhiLambda, BetaLambda = BetaLambda, LogPosterior = lpost))
    
}



# xfile function: function that estimates the model
  # arg:
  # return:
xfile = function(y, x=NULL, w=NULL, hyperparameter, 
                 data_distribution ="gaussian", 
                 Sigma = NULL, tolerance = 10^(-8), max_it = 500, max_nfact = min(ncol(y), nrow(x)),
                 seed = 1111){
  n = nrow(y)
  p = ncol(y)
  alpha_ratio = (1+hyperparameter$alpha)/hyperparameter$alpha
  Factors = list()
  if(is.null(Sigma)&!(data_distribution %in% c("binary", "Binary", "probit", "Probit"))){
    listfunction = sigmaNullFunctions(hyperparameter$a_sigma, hyperparameter$b_sigma)
  }else{
    listfunction = sigmaKnownFunctions()
  }
  
  Factors[[1]] = xfileCoordinateMM(y=y, n, p, x=x, w=w, M_h = matrix(0,n,p), hyperparameter = hyperparameter,
                          data_distribution = data_distribution, 
                          listfunction = listfunction, tolerance = tolerance, max_it = max_it,
                          seed = seed)
  n_it = length(Factors[[1]]$Vartheta)
  M_h = F_h = Factors[[1]]$F_h
  Zresidual = Factors[[1]]$Z[n_it,,]-F_h
  Vartheta = Factors[[1]]$Vartheta[n_it]
  Eta = Factors[[1]]$Eta[n_it,]
  PhiEta = Factors[[1]]$PhiEta[n_it,]
  BetaEta = Factors[[1]]$BetaEta[n_it,]
  Lambda = Factors[[1]]$Lambda[n_it,]
  PhiLambda = Factors[[1]]$PhiLambda[n_it,]
  BetaLambda = Factors[[1]]$BetaLambda[n_it,]
  LogPosterior = c(Factors[[1]]$LogPosterior[1], Factors[[1]]$LogPosterior[n_it])
  
  cat("Factor 1 estimated \n")
  
  h=2
  end=F
  while((end==F)&(h<(max_nfact+1))){
    
    Factors[[h]] = xfileCoordinateMM(y=y, x=x, w=w, M_h = M_h, hyperparameter = hyperparameter,
                      data_distribution = data_distribution, 
                      listfunction = listfunction, tolerance = tolerance, max_it = max_it,
                      seed = seed)
    n_it = length(Factors[[h]]$Vartheta)
    
    l0 = listfunction$loglikelihood(Zresidual, mu=0)
    Zresidualtmp = Factors[[h]]$Z[n_it,,]
    lF = listfunction$loglikelihood(Zresidualtmp, mu=Factors[[h]]$F_h)

    if(log( alpha_ratio^(h+1) -1) < sum(lF-l0, na.rm = T)){
      M_h = M_h + Factors[[h]]$F_h
      Zresidual = Zresidualtmp-Factors[[h]]$F_h
      Vartheta = c(Vartheta, Factors[[h]]$Vartheta[n_it])
      Eta = cbind(Eta, Factors[[h]]$Eta[n_it,])
      PhiEta = cbind(PhiEta, Factors[[h]]$PhiEta[n_it,])
      BetaEta = cbind(BetaEta, Factors[[h]]$BetaEta[n_it,])
      Lambda = cbind(Lambda, Factors[[h]]$Lambda[n_it,]) 
      PhiLambda = cbind(PhiLambda, Factors[[h]]$PhiLambda[n_it,])
      BetaLambda = cbind(BetaLambda, Factors[[h]]$BetaLambda[n_it,])
      LogPosterior = c(LogPosterior, Factors[[h]]$LogPosterior[n_it])
      
      cat("Factor ",h," estimated \n")
      h = h+1
    }else{end=T}
    
  }
  
  # names
  rownames(Eta) = rownames(PhiEta) = rownames(y)
  rownames(Lambda) = rownames(PhiLambda) = colnames(y)
  rownames(BetaEta) = colnames(x)
  rownames(BetaLambda) = colnames(w)
  
  return(list(Model = M_h, LogPosterior = LogPosterior, k = h-1, Vartheta = Vartheta, GaussianResidual = Zresidual,
         Eta = Eta, PhiEta = PhiEta, BetaEta = BetaEta,
         Lambda = Lambda, PhiLambda = PhiLambda, BetaLambda = BetaLambda))
}



# xfile simulation
#
#
xfileSimulation = function(n, p, k, c_n, c_p, epsilon = 0.1, 
                           x = matrix(1,n,1), w = matrix(1, p,1),
                           data_distribution = "gaussian",
                           data_generator_process = "additive",
                           NA_frac = 0.25){
  qx = ncol(x)
  qw = ncol(w)
  betaeta = matrix(rnorm(qx*k,0,1), qx, k)
  betaeta = betaeta+ 2*epsilon*sign(betaeta)
  betalambda = matrix(rnorm(qw*k,0,1), qw,k)
  betalambda = betalambda+ 2*epsilon*sign(betalambda)
  
  eta = matrix(rnorm(n*k),n,k)
  eta = eta+ epsilon*2*sign(eta)
  if(data_generator_process == "additive"){
    Eta0 = eta+(pmax(drop(x%*%betaeta),0)+epsilon)
  }else if(data_generator_process == "multiplicative"){
    Eta0 = eta*(pmax(drop(x%*%betaeta),0)+epsilon)
  }else{
    stop("The value of data_generator_process must be either multiplicative or additive")  
  }
  for(h in 1:k){
    Eta0[order(abs(Eta0[,h]))[1:round((1-c_n)*n)], h] = 0
  }
  Eta0 = Eta0[,order(apply(Eta0,2,var) , decreasing=T)] 
  
  lambda = matrix(rnorm(p*k),p,k)
  lambda = lambda+ epsilon*2*sign(lambda)
  if(data_generator_process == "additive"){
    Lambda0 = lambda+(pmax(drop(w%*%betalambda),0)+epsilon)
  }else if(data_generator_process == "multiplicative"){
    Lambda0 = lambda*(pmax(drop(w%*%betalambda),0)+epsilon)
  }
  for(h in 1:k){
    Lambda0[order(abs(Lambda0[,h]))[1:round((1-c_p)*p)],h] = 0
  }
  Lambda0 = Lambda0[,order(apply(Lambda0,2,var) , decreasing=T)] 
  
  z = tcrossprod(Eta0,Lambda0)+ matrix(rnorm(n*p, 0, sd=0.25),n,p)
  
  if(data_distribution %in% c("gaussian", "Gaussian", "normal", "Normal")){
    y = z
  }else if(data_distribution %in% c("binary", "Binary", "probit", "Probit")){
    y = round(z/max(z)+0.5)
  }else if(data_distribution %in% c("censored", "Censored", "Truncated", "Truncated")){
    y = pmax(z,0)
  }
  
  y_na = y
  y_na[ sample(n*p, round(n*p*NA_frac))] = NA
  na_col = colSums(is.na(y_na))
  na_row = rowSums(is.na(y_na))
  wnc = which(na_col == n)
  if(length(wnc)>0){
    swnc = sample(1:n,wnc)
    y_na[swnc, wnc] = y[swnc, wnc]
  }
  wnr = which(na_col == p)
  if(length(wnr)>0){
    swnr = sample(1:p,wnr)
    y_na[swnr, wnr] = y[swnr, wnr]
  }

  return( list(y = y, y_na = y_na, z0 = z, Eta0 = Eta0, Lambda0 = Lambda0,
               BetaEta0 = betaeta, BetaLambda0 = betalambda))
}



sigmaKnownFunctions=function(){
  
  WetainvUpdate = function()
    DetaUpdate = function(){}
    DbetaetaUpdate = function(){}
    WbetaetaUpdate = function(){}
    constant = 
      
      Dlambdaupdate = function(){}
    Wlambdaupdate = function(){}
    DbetalambdaUpdate = function(){}
    WbetalambdaUpdate = function(){}
    
    return(list(DetaUpdate, DlambdaUpdate, DbetaetaUpdate, WbetaetaUpdate,
                WetaUpdate, WlambdaUpdate, DbetalambdaUpdate, WbetalambdaUpdate, constant
    ))
}
