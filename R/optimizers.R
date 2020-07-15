# Algorithms to optimize the GLM-PCA objective function.

#Global variables
NB_THETA_MAX<-1e4 
#negative binomial theta values larger than this are truncated to this value
#motivation: large nb_theta makes it essentially poisson, no point in estimating if beyond this point.

check_divergence<-function(curr,alg=c("avagrad","avagrad_stochastic","fisher"),
                           ctl_param){
  #ctl is the value of 'lr' (if avagrad) or 'penalty' (if fisher)
  if(!is.finite(curr)){
    alg<-match.arg(alg)
    msg1<-"Numerical divergence (deviance no longer finite), try"
    msg3<-"control parameter to improve stability of optimization."
    msg4<-paste("Current control parameter value:",ctl_param)
    if(alg=="avagrad"){
      msg2<-"decreasing the learning rate (lr)"
    } else if(alg=="avagrad_stochastic"){
      msg2<-"increasing the minibatch size or decreasing the learning rate (lr)"
    } else { #alg=="fisher"
      msg2<-"increasing the penalty"
    }
    stop_custom("error_glmpca_divergence",paste(msg1,msg2,msg3,msg4))
  }
}

#' @importFrom utils tail
check_dev_decr<-function(dev){
  if(tail(dev,1)>dev[1]){
    msg1<-"Poor model fit (final deviance higher than initial deviance)."
    msg2<-"Try modifying control parameters to improve optimization."
    stop_custom("error_glmpca_dev_incr",paste(msg1,msg2))
  }
  dev
}

# check_convergence_smoothed<-function(x,tol){
#   #x is a numeric vector of length >=2
#   #tries to average over recent noisy objective function values
#   #to determine convergence
#   #if length(x)==2, mad(x,constant=1) equiv to abs(x[2]-x[1])/2
#   mad(x,constant=1)/(0.1+abs(median(x))) < tol
# }

check_convergence<-function(x1,x2,tol,thresh=NULL){
  #x1 is old objective function value
  #x2 is new objective function value
  below_thresh<- !is.null(thresh) && x2<thresh #usually is FALSE
  below_thresh || abs(x2-x1)/(0.1+abs(x1)) < tol
}

print_status<-function(curr,iter,nb_theta=NULL){
  dev_format<-format(curr,scientific=TRUE,digits=4)
  msg<-paste0("Iteration: ",iter," | deviance=",dev_format)
  if(length(nb_theta)==1){ msg<-paste0(msg," | nb_theta: ",signif(nb_theta,3)) }
  message(msg) 
}

nb_theta_infograd<-function(Y,Mu,th,grad_only=TRUE){
  #given count data y and predicted means mu>0, and a vector of negative binomial theta "th" params
  #note: Y, Mu are assumed to be matrices whose rows are features
  #return the gradient and observed information of the log likelihood with respect to the negative binomial theta (th) params
  #note this uses observed rather than expected information
  #dtheta/du=e^u=theta
  #d2theta/du2=theta
  if(length(th)==1){
    fam<-"nb"
    reduce_func<-sum
  } else {
    fam<-"nb2"
    stopifnot(length(th)==nrow(Y))
    reduce_func<-rowSums
  }
  #dL/dtheta*dtheta/du
  #note: th+Y uses recycling! We assume length(th)==nrow(Y)
  ig<-list(grad=th*reduce_func(digamma(th+Y)-digamma(th)+log(th)+1-log(th+Mu)-(Y+th)/(Mu+th)))
  if(grad_only){ 
    return(ig) 
  } else {
    #d^2L/dtheta^2 * (dtheta/du)^2
    info1<- -th^2*reduce_func(trigamma(th+Mu)-trigamma(th)+1/th-2/(Mu+th)+(Y+th)/(Mu+th)^2)  
    #dL/dtheta*d^2theta/du^2 = grad
    ig$info<- info1-ig$grad
    return(ig)
  }
}

est_nb_theta_fisher<-function(Y,Mu,th,logscale_prior=c(0,1)){
  #use Newton's Method to update theta vector based on the negative binomial likelihood
  #logscale_prior is a lognormal prior for th for regularization.
  #let u=log(theta). We use the prior u~N(u0,1/t0) as penalty
  #equivalently we assume theta~lognormal(u0,1/t0) so the mode is at exp(mu)
  #to remove the effect of the prior set logscale_prior=c(0,0)
  #note the t0 here is the inverse of the Gaussian variance parameter (ie the precision)
  u<-log(th)
  ig<-nb_theta_infograd(Y,Mu,th,grad_only=FALSE)
  u0<-logscale_prior[1]; t0<-logscale_prior[2]
  grad<-ig$grad - t0*(u-u0)
  info<-ig$info + t0
  #L2 penalty on u=log(th)
  exp(u+grad/info)
}

fisher_optimizer<-function(Y,U,V,uid,vid,ctl,gf,rfunc,offsets){
  #Y: the data matrix
  #U: initialized factors matrix, including all column covariates & coefficients
  #V: initialized loadings matrix, including all row covariates & coefficients
  #uid: which columns of U contain updateable parameters?
  #vid: which columns of V contain updateable parameters?
  #ctl: list of optimization control params: penalty,maxIter,minIter,and tol
  #gf: object of type glmpca_family
  #...for fam 'poi','nb','nb2' gf contains the 'offsets' 
  #...for fam 'binom' offsets evaluates to NULL, can still be passed to rfunc
  #rfunc: a function that computes the linear predictor ...
  #...of the GLM from U,V, and offsets.
  
  #lid: which columns of U,V contain the unsupervised factors & loadings
  lid<-intersect(uid,vid)
  stopifnot(all(c("penalty","maxIter","minIter","tol") %in% names(ctl)))
  sz<-if(gf$glmpca_fam=="binom"){ gf$binom_n } else { NULL }
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){
    #rmse[t]<-sd(Y-ilfunc(rfunc(U,V,offsets)))
    dev[t]<-gf$dev_func(Y,rfunc(U,V,offsets),sz=sz)
    check_divergence(dev[t],"fisher",ctl$penalty)
    if(ctl$verbose){print_status(dev[t],t,gf$nb_theta)}
    if(t>ctl$minIter && check_convergence(dev[t-1],dev[t],ctl$tol,ctl$minDev)){
      break
    }
    # fisher scoring optimizer, slow b/c recomputes gradient for each dimension.
    # (k %in% lid) ensures no penalty on regression coefficients:
    # note: ig$infograd gives gradient for the log-likelihood
    # the grad for the deviance has opposite sign (-grads).
    # note this does not affect diagonal Fisher scoring since the signs cancel.
    for(k in uid){
      ig<- gf$infograd(Y, rfunc(U,V,offsets), sz=sz, grad_only=FALSE)
      grads<- crossprod(ig$grad, V[,k]) - ctl$penalty*U[,k]*(k %in% lid)
      infos<- crossprod(ig$info, V[,k]^2) + ctl$penalty*(k %in% lid)
      U[,k]<-U[,k]+grads/infos
    }
    for(k in vid){
      ig<- gf$infograd(Y, rfunc(U,V,offsets), sz=sz, grad_only=FALSE)
      grads<- (ig$grad)%*%U[,k] - ctl$penalty*V[,k]*(k %in% lid)
      infos<- (ig$info) %*% U[,k]^2 + ctl$penalty*(k %in% lid)
      V[,k]<-V[,k]+grads/infos
    }
    
    # alternative: fisher scoring on full blocks, very unstable numerically
    # ig<- gf$infograd(Y,rfunc(U,V,offsets))
    # grad_u<- crossprod(ig$grad, V[,uid]) - penalty*U[,uid]*(uid %in% lid)
    # info_u<- crossprod(ig$info, V[,uid]^2) + penalty*(uid %in% lid)
    # grad_v<- (ig$grad)%*%U[,vid] - penalty*V[,vid]*(vid %in% lid)
    # info_v<- (ig$info) %*% U[,vid]^2 + penalty*(vid %in% lid)
    # U[,uid]<-U[,uid]+lr*grad_u/info_u
    # V[,vid]<-V[,vid]+lr*grad_v/info_v
    
    if(gf$glmpca_fam %in% c("nb","nb2")){
      #pmin here is to prevent extremely large values, which don't make much difference in the model fitting.
      nb_theta<-est_nb_theta_fisher(Y, gf$linkinv(rfunc(U,V,offsets)), gf$nb_theta)
      gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,nb_theta))
    }
  }
  list(U=U, V=V, dev=check_dev_decr(dev[1:t]), gf=gf)
}

#in the avagrad paper, lr=alpha (step size)
#eps=epsilon controls the "adaptivity"
#betas=beta1,beta2 are inherited from the Adam method
#usually only need to tune lr for avagrad.
#avagrad paper: https://arxiv.org/abs/1912.01823
#avagrad pytorch implementation: ...
#...https://github.com/lolemacs/avagrad/blob/master/optimizers.py
#avagrad illustrated on toy function ...
#https://github.com/willtownes/optimization-practice/blob/master/practice.Rmd#L159

avagrad_optimizer<-function(Y,U,V,uid,vid,ctl,gf,rfunc,offsets){
  #Y: the data matrix
  #U: initialized factors matrix, including all column covariates & coefficients
  #V: initialized loadings matrix, including all row covariates & coefficients
  #uid: which columns of U contain updateable parameters?
  #vid: which columns of V contain updateable parameters?
  #ctl: list of optimization control params: maxIter, minIter, tol, ...
  #... ctl$betas, ctl$epsilon, ctl$lr are Adam tuning params
  #gf: object of type glmpca_family
  #rfunc: a function that computes the linear predictor ...
  #...of the GLM from U,V, and offsets.
  
  #lid: which columns of U,V contain the unsupervised factors & loadings
  lid<-intersect(uid,vid)
  #avagrad initialization: prefix m_ refers to momentum, v_ refers to variance
  m_u<-v_u<-matrix(0,nrow=ncol(Y),ncol=length(uid))
  m_v<-v_v<-matrix(0,nrow=nrow(Y),ncol=length(vid))
  k<-length(m_u)+length(m_v) #total number of parameters being updated
  # if(gf$glmpca_fam %in% c("nb","nb2")){
  #   m_th<-v_th<- 0*gf$nb_theta
  #   k<- k+length(gf$nb_theta)
  # }
  sz<-if(gf$glmpca_fam=="binom"){ gf$binom_n } else { NULL }
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){
    dev[t]<-gf$dev_func(Y, rfunc(U,V,offsets), sz=sz)
    check_divergence(dev[t],"avagrad",ctl$lr)
    if(ctl$verbose){print_status(dev[t],t,gf$nb_theta)}
    if(t>ctl$minIter && check_convergence(dev[t-1],dev[t],ctl$tol,ctl$minDev)){
      break
    }
    # avagrad optimizer, does not need ig$info, sensitive to learning rate, ...
    #... no penalty needed but has additional tuning params betas, eps, lr
    ig<-gf$infograd(Y, rfunc(U,V,offsets), sz=sz, grad_only=TRUE)
    #ig$grad contains gradient of log-likelihood. Deviance has opposite sign.
    g_u<- -crossprod(ig$grad, V[,uid]) #+ ctl$penalty*t(t(U[,uid])*(uid %in% lid))
    g_v<- -(ig$grad)%*%U[,vid] #+ ctl$penalty*t(t(V[,vid])*(vid %in% lid))
    m_u<-ctl$betas[1]*m_u+(1-ctl$betas[1])*g_u
    m_v<-ctl$betas[1]*m_v+(1-ctl$betas[1])*g_v
    eta_u<-1/(ctl$epsilon+sqrt(v_u))
    eta_v<-1/(ctl$epsilon+sqrt(v_v))
    U[,uid]<-U[,uid] - ctl$lr*(eta_u/l2norm(eta_u/sqrt(k)))*m_u
    V[,vid]<-V[,vid] - ctl$lr*(eta_v/l2norm(eta_v/sqrt(k)))*m_v
    v_u<-ctl$betas[2]*v_u+(1-ctl$betas[2])*g_u^2
    v_v<-ctl$betas[2]*v_v+(1-ctl$betas[2])*g_v^2
    
    if(gf$glmpca_fam %in% c("nb","nb2")){  
      #note: all the g,m,eta, terms here are with respect to log(nb_theta)
      th<-gf$nb_theta
      # lth<-log(th)
      # ig_th<- nb_theta_infograd(Y, gf$linkinv(rfunc(U,V,offsets)), th, grad_only=TRUE)
      # g_th<- -(ig_th$grad)
      # m_th<- ctl$betas[1]*m_th + (1-ctl$betas[1])*g_th
      # eta_th<- 1/(ctl$epsilon+sqrt(v_th))
      # lth<- lth - ctl$lr*(eta_th/l2norm(eta_th/sqrt(k)))*m_th
      # v_th<- ctl$betas[2]*v_th+(1-ctl$betas[2])*g_th^2
      # gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,exp(lth)))
      nb_theta<-est_nb_theta_fisher(Y, gf$linkinv(rfunc(U,V,offsets)), th)
      gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,nb_theta))
    }
  }
  list(U=U, V=V, dev=check_dev_decr(dev[1:t]), gf=gf)
}

#this one updates loadings V after each epoch, equivalent to full gradient
avagrad_memoized_optimizer<-function(Y,U,V,uid,vid,ctl,gf,rfunc,offsets){
  #Y: the data matrix
  #U: initialized factors matrix, including all column covariates & coefficients
  #V: initialized loadings matrix, including all row covariates & coefficients
  #uid: which columns of U contain updateable parameters?
  #vid: which columns of V contain updateable parameters?
  #ctl: list of optimization control params: maxIter, minIter, tol, ...
  #... ctl$betas, ctl$epsilon, ctl$lr are Adam tuning params
  #gf: object of type glmpca_family
  #rfunc: a function that computes the linear predictor ...
  #...of the GLM from U,V, and offsets.
  N<-ncol(Y); J<-nrow(Y); Kv<-length(vid)
  #lid: which columns of U,V contain the unsupervised factors & loadings
  lid<-intersect(uid,vid)
  #avagrad initialization: prefix m_ refers to momentum, v_ refers to variance
  m_u<-v_u<-matrix(0,nrow=N,ncol=length(uid))
  m_v<-v_v<-matrix(0,nrow=J,ncol=Kv)
  k<-length(m_u)+length(m_v) #total number of parameters being updated
  # if(gf$glmpca_fam %in% c("nb","nb2")){ m_th<-v_th<- 0*gf$nb_theta }
  
  #initialize minibatches (fixed for all epochs) and sufficient statistics
  mb_list<-create_minibatches(N,ctl$batch_size,randomize=TRUE)
  B<-length(mb_list)
  #3D array, grad of V will be rowSums(ss,dims=2) (ie sum on the last dim)
  S<-array(0,dim=c(J,Kv,B)) 
  #g_v<-m_v
  sz<-if(gf$glmpca_fam=="binom"){ gf$binom_n } else { NULL }
  
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){ #one iteration = 1 epoch
    devb<-0 #accumulator for deviance across each minibatch
    for(b in 1:B){
      mb<-mb_list[[b]]
      Ymb<-as.matrix(Y[,mb]) #assume this can fit in memory.
      szb<-sz[mb]
      
      # avagrad optimizer, does not need ig$info, sensitive to learning rate, ...
      #... no penalty needed but has additional tuning params betas, eps, lr
      ig<-gf$infograd(Ymb, rfunc(U[mb,],V,offsets[mb]), sz=szb, grad_only=TRUE)
      #ig$grad contains gradient of log-likelihood. Deviance has opposite sign.
      g_u<- -crossprod(ig$grad, V[,uid])
      S[,,b]<-(-(ig$grad)%*%U[mb,vid])
      m_u[mb,]<-ctl$betas[1]*m_u[mb,]+(1-ctl$betas[1])*g_u
      eta_u<-1/(ctl$epsilon+sqrt(v_u[mb,]))
      U[mb,uid]<-U[mb,uid] - ctl$lr*(eta_u/l2norm(eta_u/sqrt(k)))*m_u[mb,]
      v_u[mb,]<-ctl$betas[2]*v_u[mb,]+(1-ctl$betas[2])*g_u^2
      
      R<-rfunc(U[mb,],V,offsets[mb])
      devb<-devb + gf$dev_func(Ymb,R,sz=szb)
      check_divergence(devb,"avagrad",ctl$lr)
      if(gf$glmpca_fam == "nb"){
        th<-gf$nb_theta
        nb_theta<-est_nb_theta_fisher(Ymb, gf$linkinv(R), th)
        #take weighted average of previous and current theta estimate
        #weight is set so that equivalent of a full fisher update once per epoch
        wt<-1/B
        nb_theta<- exp(wt*log(nb_theta)+(1-wt)*log(th))
        gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,nb_theta))
      } #end of nb theta update
    } #end of loop over minibatches (end of epoch)
    
    #update global params (loadings) V at end of epoch
    g_v<-rowSums(S,dims=2) #sum over 3D array accumulating sufficient stats
    m_v<-ctl$betas[1]*m_v+(1-ctl$betas[1])*g_v
    eta_v<-1/(ctl$epsilon+sqrt(v_v))
    V[,vid]<-V[,vid] - ctl$lr*(eta_v/l2norm(eta_v/sqrt(k)))*m_v
    v_v<-ctl$betas[2]*v_v+(1-ctl$betas[2])*g_v^2
    
    #assess convergence after end of each epoch
    dev[t]<-devb
    if(ctl$verbose){print_status(dev[t],t,gf$nb_theta)}
    if(t>ctl$minIter && check_convergence(dev[t-1],dev[t],ctl$tol,ctl$minDev)){
      break
    }
  }
  list(U=U, V=V, dev=check_dev_decr(dev[1:t]), gf=gf)
}

#this one updates loadings V after each minibatch (ie more often)
avagrad_memoized_optimizer2<-function(Y,U,V,uid,vid,ctl,gf,rfunc,offsets){
  #Y: the data matrix
  #U: initialized factors matrix, including all column covariates & coefficients
  #V: initialized loadings matrix, including all row covariates & coefficients
  #uid: which columns of U contain updateable parameters?
  #vid: which columns of V contain updateable parameters?
  #ctl: list of optimization control params: maxIter, minIter, tol, ...
  #... ctl$betas, ctl$epsilon, ctl$lr are Adam tuning params
  #gf: object of type glmpca_family
  #rfunc: a function that computes the linear predictor ...
  #...of the GLM from U,V, and offsets.
  N<-ncol(Y); J<-nrow(Y); Kv<-length(vid)
  #lid: which columns of U,V contain the unsupervised factors & loadings
  lid<-intersect(uid,vid)
  #avagrad initialization: prefix m_ refers to momentum, v_ refers to variance
  m_u<-v_u<-matrix(0,nrow=N,ncol=length(uid))
  m_v<-v_v<-matrix(0,nrow=J,ncol=Kv)
  # if(gf$glmpca_fam %in% c("nb","nb2")){ m_th<-v_th<- 0*gf$nb_theta }
  
  #initialize minibatches (fixed for all epochs) and sufficient statistics
  mb_list<-create_minibatches(N,ctl$batch_size,randomize=TRUE)
  B<-length(mb_list)
  #3D array, grad of V will be rowSums(S,dims=2) (ie sum on the last dim)
  S<-array(0,dim=c(J,Kv,B)) 
  g_v<-m_v
  sz<-if(gf$glmpca_fam=="binom"){ gf$binom_n } else { NULL }
  
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){ #one iteration = 1 epoch
    devb<-0 #accumulator for deviance across each minibatch
    for(b in 1:B){
      mb<-mb_list[[b]]
      mb_len<-length(mb)
      k<-length(m_v) + mb_len*ncol(m_u) #total number of parameters being updated
      # if(gf$glmpca_fam %in% c("nb","nb2")){
      #   k<- k+length(gf$nb_theta)
      # }
      #adj_factor<-N/mb_len
      Ymb<-as.matrix(Y[,mb]) #assume this can fit in memory.
      szb<-sz[mb]
      
      # avagrad optimizer, does not need ig$info, sensitive to learning rate, ...
      #... no penalty needed but has additional tuning params betas, eps, lr
      ig<-gf$infograd(Ymb, rfunc(U[mb,],V,offsets[mb]), sz=szb, grad_only=TRUE)
      #ig$grad contains gradient of log-likelihood. Deviance has opposite sign.
      g_u<- -crossprod(ig$grad, V[,uid])
      if(t>1){ g_v<- g_v - S[,,b] } #remove suff stat from previous epoch
      S[,,b]<-(-(ig$grad)%*%U[mb,vid])
      if(t>1){ #full gradient already computed
        g_v<- g_v + S[,,b]
      } else { #first epoch, rescale cumulative minibatches as if stoch grad
        adj_factor<-N/length(unlist(mb_list[1:b]))
        g_v<- adj_factor*rowSums(S[,,1:b,drop=FALSE],dims=2)
      }
      m_u[mb,]<-ctl$betas[1]*m_u[mb,]+(1-ctl$betas[1])*g_u
      m_v<-ctl$betas[1]*m_v+(1-ctl$betas[1])*g_v
      eta_u<-1/(ctl$epsilon+sqrt(v_u[mb,]))
      eta_v<-1/(ctl$epsilon+sqrt(v_v))
      U[mb,uid]<-U[mb,uid] - ctl$lr*(eta_u/l2norm(eta_u/sqrt(k)))*m_u[mb,]
      V[,vid]<-V[,vid] - ctl$lr*(eta_v/l2norm(eta_v/sqrt(k)))*m_v
      v_u[mb,]<-ctl$betas[2]*v_u[mb,]+(1-ctl$betas[2])*g_u^2
      v_v<-ctl$betas[2]*v_v+(1-ctl$betas[2])*g_v^2

      R<-rfunc(U[mb,],V,offsets[mb])
      devb<-devb + gf$dev_func(Ymb,R,sz=szb)
      check_divergence(devb,"avagrad",ctl$lr)
      if(gf$glmpca_fam == "nb"){
        th<-gf$nb_theta
        nb_theta<-est_nb_theta_fisher(Ymb, gf$linkinv(R), th)
        #take weighted average of previous and current theta estimate
        #weight is set so that equivalent of a full fisher update once per epoch
        wt<-1/B
        nb_theta<- exp(wt*log(nb_theta)+(1-wt)*log(th))
        gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,nb_theta))
      } #end of nb theta update
    } #end of loop over minibatches (end of epoch)
    
    #assess convergence after end of each epoch
    dev[t]<-devb
    if(ctl$verbose){print_status(dev[t],t,gf$nb_theta)}
    if(t>ctl$minIter && check_convergence(dev[t-1],dev[t],ctl$tol,ctl$minDev)){
      break
    }
  }
  list(U=U, V=V, dev=check_dev_decr(dev[1:t]), gf=gf)
}

#' @importFrom stats fitted lm
#' @importFrom utils tail
avagrad_stochastic_optimizer<-function(Y,U,V,uid,vid,ctl,gf,rfunc,offsets){
  #Y: the data matrix
  #U: initialized factors matrix, including all column covariates & coefficients
  #V: initialized loadings matrix, including all row covariates & coefficients
  #uid: which columns of U contain updateable parameters?
  #vid: which columns of V contain updateable parameters?
  #ctl: list of optimization control params: maxIter, minIter, tol, ...
  #... ctl$betas, ctl$epsilon, ctl$lr are Adam tuning params
  #gf: object of type glmpca_family
  #rfunc: a function that computes the linear predictor ...
  #...of the GLM from U,V, and offsets.
  stopifnot(ctl$minIter>10)
  N<-ncol(Y)
  #lid: which columns of U,V contain the unsupervised factors & loadings
  lid<-intersect(uid,vid)
  #avagrad initialization: prefix m_ refers to momentum, v_ refers to variance
  m_u<-v_u<-matrix(0,nrow=N,ncol=length(uid))
  m_v<-v_v<-matrix(0,nrow=nrow(Y),ncol=length(vid))
  # if(gf$glmpca_fam %in% c("nb","nb2")){ m_th<-v_th<- 0*gf$nb_theta }
  sz<-if(gf$glmpca_fam=="binom"){ gf$binom_n } else { NULL }
  
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  dev_smooth<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){ #one iteration = 1 epoch
    #create minibatch indices
    mb_list<-create_minibatches(N,ctl$batch_size,randomize=TRUE)
    B<-length(mb_list)
    #b<- ((t-1) %% B) + 1 #1,2,3,4,...,B,1,2,3,4,...,B,....
    for(b in 1:B){
      mb<-mb_list[[b]]
      mb_len<-length(mb)
      k<-length(m_v) + mb_len*ncol(m_u) #total number of parameters being updated
      # if(gf$glmpca_fam %in% c("nb","nb2")){
      #   k<- k+length(gf$nb_theta)
      # }
      adj_factor<-N/mb_len
      Ymb<-as.matrix(Y[,mb]) #assume this can fit in memory.
      szb<-sz[mb]
      
      # avagrad optimizer, does not need ig$info, sensitive to learning rate, ...
      #... no penalty needed but has additional tuning params betas, eps, lr
      ig<-gf$infograd(Ymb, rfunc(U[mb,],V,offsets[mb]), sz=szb, grad_only=TRUE)
      #ig$grad contains gradient of log-likelihood. Deviance has opposite sign.
      g_u<- -crossprod(ig$grad, V[,uid]) #+ ctl$penalty*t(t(U[mb,uid])*(uid %in% lid))
      g_v<- adj_factor*(-(ig$grad)%*%U[mb,vid]) #+ ctl$penalty*t(t(V[,vid])*(vid %in% lid))
      m_u[mb,]<-ctl$betas[1]*m_u[mb,]+(1-ctl$betas[1])*g_u
      m_v<-ctl$betas[1]*m_v+(1-ctl$betas[1])*g_v
      eta_u<-1/(ctl$epsilon+sqrt(v_u[mb,]))
      eta_v<-1/(ctl$epsilon+sqrt(v_v))
      U[mb,uid]<-U[mb,uid] - ctl$lr*(eta_u/l2norm(eta_u/sqrt(k)))*m_u[mb,]
      V[,vid]<-V[,vid] - ctl$lr*(eta_v/l2norm(eta_v/sqrt(k)))*m_v
      v_u[mb,]<-ctl$betas[2]*v_u[mb,]+(1-ctl$betas[2])*g_u^2
      v_v<-ctl$betas[2]*v_v+(1-ctl$betas[2])*g_v^2
      
      R<-rfunc(U[mb,],V,offsets[mb])
      # if(gf$glmpca_fam %in% c("nb","nb2")){
      if(gf$glmpca_fam == "nb"){
        th<-gf$nb_theta
        #note: all the g,m,eta, terms here are with respect to log(nb_theta)
        # lth<-log(th)
        # ig_th<- nb_theta_infograd(Ymb, gf$linkinv(rfunc(U[mb,],V,offsets[mb])), th, grad_only=TRUE)
        # g_th<- adj_factor*(-ig_th$grad)
        # m_th<- ctl$betas[1]*m_th + (1-ctl$betas[1])*g_th
        # eta_th<- 1/(ctl$epsilon+sqrt(v_th))
        # lth<- lth - ctl$lr*(eta_th/l2norm(eta_th/sqrt(k)))*m_th
        # v_th<- ctl$betas[2]*v_th+(1-ctl$betas[2])*g_th^2
        # gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,exp(lth)))
        nb_theta<-est_nb_theta_fisher(Ymb, gf$linkinv(R), th)
        #take weighted average of previous and current theta estimate
        #weight is set so that equivalent of a full fisher update once per epoch
        wt<-1/B
        nb_theta<- exp(wt*log(nb_theta)+(1-wt)*log(th))
        gf<-glmpca_family(gf$glmpca_fam, nb_theta=pmin(NB_THETA_MAX,nb_theta))
      } #end of nb theta update
    } #end of loop over minibatches (end of epoch)
    
    #assess convergence after end of each epoch
    dev[t]<-adj_factor*gf$dev_func(Ymb,R,sz=szb)
    check_divergence(dev[t],"avagrad",ctl$lr)
    j<-seq.int(max(1,t-10+1),t)
    #dev_smooth[t]<-exp(mean(log(dev[j]))) #mean(dev[j])
    dev_smooth[t]<-exp(tail(fitted(lm(log(dev[j])~j)),1))
    # if(t>B){
    #   fit<-StructTS(log(dev[1:t]),type="level")
    #   dev_smooth[t]<-exp(tail(fitted(fit),1))
    # } else {
    #   j<-seq.int(max(1,t-B+1),t)
    #   dev_smooth[t]<-exp(tail(fitted(lm(log(dev[j])~j)),1))
    # }
    if(ctl$verbose){print_status(dev[t],t,gf$nb_theta)}
    if(t>ctl$minIter && (dev_smooth[t]< dev[1])){
      #ensure at least one full pass through the data
      if(check_convergence(dev_smooth[t-1],dev_smooth[t],ctl$tol,ctl$minDev)){
       break
      }
      # if(fit$coef["level"]<ctl$tol){ break }
    }
  }
  list(U=U, V=V, dev=check_dev_decr(dev[1:t]), gf=gf, 
       dev_smooth=check_dev_decr(dev_smooth[1:t]))
}
