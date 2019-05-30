# GLM-PCA
#source("./algs/ortho.R")

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  # apply(x,2,norm)
  sqrt(colSums(x^2))
}

ortho<-function(U,V,A,X=rep(1,nrow(U)),G=0,Z=0,ret=c("m","df")){
  #U is NxL matrix of cell factors, V is JxL matrix of loadings onto genes
  #X is NxKo matrix of cell specific covariates
  #A is JxKo matrix of coefficients of X
  #Z is JxKf matrix of gene specific covariates
  #G is NxKf matrix of coefficients of Z
  #ret= return either data frame or matrix format
  #imputed expression: E[Y] = g^{-1}(R) where R = AX'+ZG'+VU'
  ret<-match.arg(ret)
  L<-ncol(U)
  X<-as.matrix(X,nrow=nrow(U))
  if(is.null(Z)){ Z<-0 }
  Z<-as.matrix(Z,nrow=nrow(V))
  if(all(G==0)){ G<-NULL }
  #we assume A is not null or zero
  #remove correlation between U and A
  #at minimum, this will cause factors to have mean zero
  reg<-lm(U~X-1)
  factors<-residuals(reg)
  A<-A+tcrossprod(V,coef(reg))
  #remove correlation between V and G
  if(is.null(G)){
    loadings<-V
  } else { #G is not empty
    reg<-lm(V~Z-1)
    loadings<-residuals(reg)
    G<-G+tcrossprod(factors,coef(reg))
  }
  #rotate factors to make loadings orthornormal
  svdres<-svd(loadings)
  loadings<-svdres$u
  factors<-t(t(factors%*%svdres$v)*svdres$d)
  #arrange latent dimensions in decreasing L2 norm
  o<-order(colNorms(factors),decreasing=TRUE)
  factors<-factors[,o]
  loadings<-loadings[,o]
  colnames(loadings)<-colnames(factors)<-paste0("dim",1:L)
  if(ret=="df"){
    loadings<-as.data.frame(loadings)
    factors<-as.data.frame(factors)
    A<-as.data.frame(A)
    if(!is.null(G)) G<-as.data.frame(G)
  }
  list(factors=factors,loadings=loadings,coefX=A,coefZ=G)
}

mat_binom_dev<-function(X,P,n){
  #binomial deviance for two matrices
  #X,P are JxN matrices
  #n is vector of length N (same as cols of X,P)
  #nz<-X>0
  X<-t(X); P<-t(P)
  term1<-sum(X*log(X/(n*P)), na.rm=TRUE)
  #nn<-x<n
  nx<-n-X
  term2<-sum(nx*log(nx/(n*(1-P))), na.rm=TRUE)
  2*(term1+term2)
}

glmpca_family<-function(fam,nb_theta=NULL,mult_n=NULL){
  #create GLM family object
  #if no offset, set offset=0
  if(fam=="poi"){
    family<-poisson()
  } else if(fam=="nb"){
    if(is.null(nb_theta)){ 
      stop("Negative binomial dispersion parameter 'nb_theta' must be specified") 
    }
    family<-MASS::negative.binomial(theta=nb_theta)
  } else if(fam %in% c("mult","bern")){
    family<-binomial()
    if(fam=="mult" && is.null(mult_n)){
      stop("Multinomial sample size parameter vector 'mult_n' must be specified")
    }
  } else {
    stop("unrecognized family type")
  }
  #variance function, determined by GLM family
  vfunc<-family$variance
  #inverse link func, mu as a function of linear predictor R
  ilfunc<-family$linkinv
  #derivative of inverse link function, dmu/dR
  hfunc<-family$mu.eta
  gf<-as.list(family)
  gf$glmpca_fam<-fam
  if(fam=="poi"){
    gf$infograd<-function(Y,R){
      M<-ilfunc(R) #ilfunc=exp
      list(grad=(Y-M),info=M)
    }
  } else if(fam=="nb"){
    gf$infograd<-function(Y,R){
      M<-ilfunc(R) #ilfunc=exp
      W<-1/vfunc(M)
      list(grad=(Y-M)*W*M, info=W*M^2)
    }
  } else if(fam=="mult"){
    gf$infograd<-function(Y,R){
      Pt<-t(ilfunc(R)) #ilfunc=expit, Pt very small probabilities
      list(grad=Y-t(mult_n*Pt), info=t(mult_n*vfunc(Pt)))
    }
  } else if(fam=="bern"){
    gf$infograd<-function(Y,R){
      P<-ilfunc(R)
      list(grad=Y-P, info=vfunc(P))
    }
  } else { #this is not actually used but keeping for future reference
    #this is most generic formula for GLM but computationally slow
    stop("invalid fam")
    gf$infograd<-function(Y,R){
      M<-ilfunc(R)
      W<-1/vfunc(M)
      H<-hfunc(R)
      list(grad=(Y-M)*W*H, info=W*H^2)
    }
  }
  #create deviance function
  if(fam=="mult"){
    gf$dev_func<-function(Y,R){
      mat_binom_dev(Y,ilfunc(R),mult_n)
    }
  } else {
    gf$dev_func<-function(Y,R){
      #Note dev.resids function gives the square of actual residuals
      #so the sum of these is the deviance
      sum(family$dev.resids(Y,ilfunc(R),1)) 
    }
  }
  class(gf)<-c("glmpca_family","family")
  gf
}

has_intercept<-function(X){
  # not used
  if(all(X[,1]==1)){
    return(TRUE)
  }
  #if 1 vector not first column, check if it is in colspace of X
  #eg if X is a saturated dummy variable encoding this is satisfied
  y<-rep(1,nrow(X))
  #if 1 vector is in colspace of X, the residuals will all be zero
  resid<-residuals(lm(y~X-1))
  max(abs(resid))<1e-12
}

remove_intercept<-function(X){
  X<-t(t(X)-colMeans(X))
  X[, colNorms(X)>1e-12, drop=FALSE]
}

glmpca_init<-function(Y,fam,sz=NULL,nb_theta=NULL){
  #create the glmpca_family object and
  #initialize the A matrix (regression coefficients of X)
  #Y is the data
  #fam is the likelihood
  #sz optional vector of size factors, default: sz=colMeans(Y) or colSums(Y)
  #sz is ignored unless fam is 'poi' or 'nb'
  if(!is.null(sz)){ stopifnot(length(sz)==ncol(Y)) }
  mult_n<- if(fam=="mult"){ colSums(Y) } else { NULL }
  gf<-glmpca_family(fam,nb_theta,mult_n)
  if(fam %in% c("poi","nb")){
    if(is.null(sz)){ sz<-colMeans(Y) } #size factors
    offsets<-gf$linkfun(sz)
    rfunc<-function(U,V){ t(offsets+tcrossprod(U,V)) } #linear predictor
    a1<-gf$linkfun(rowSums(Y)/sum(sz))
  } else {
    rfunc<-function(U,V){ tcrossprod(V,U) }
    if(fam=="mult"){ #offsets incorporated via family object
      a1<-gf$linkfun(rowSums(Y)/sum(mult_n))
    } else { #no offsets (eg, bernoulli)
      a1<-family$linkfun(rowMeans(Y))
    }
  }
  if(any(is.infinite(a1))){
    stop("Some rows were all zero, please remove them.")
  }
  list(gf=gf,rfunc=rfunc,intercepts=a1)
}

#' @title GLM-PCA
#' @description This function implements the GLM-PCA dimensionality reduction
#'   method for high-dimensional count data.
#' @name glmpca
#' 
#' @param Y matrix of count data with features as rows and observations as columns.
#' @param L the desired number of latent dimensions (integer).
#' @param fam character describing the likelihood to use for the data 
#' (poisson, negative binomial, binomial approximation to multinomial, bernoulli).
#' @param ctl a list of control parameters for optimization.
#' @param penalty the L2 penalty for the latent factors (default = 1). Regression coefficients are not penalized.
#' @param nb_theta see \code{\link[MASS]{negative.binomial}} (nb_theta -> infty = Poisson).
#' @param X a matrix of column (observations) covariates.
#' Any column with all same values (eg. 1 for intercept) will be removed.
#' This is because we force the intercept and want to avoid collinearity.
#' @param Z a matrix of row (feature) covariates, usually not needed.
#' 
#' @details The basic model is \code{R = AX'+ZG'+VU'}, where \code{E\[Y\] = M = linkinv(R)}.
#' Regression coefficients are \code{A} and \code{G}, latent factors are \code{U} and loadings are \code{V}.
#' 
#' @return a list containing:
#' 
#' @export
glmpca<-function(Y, L, fam=c("poi","nb","mult","bern"),
				 ctl = list(maxIter=1000, eps=1e-4),
				 penalty = 1, verbose = FALSE,
				 init = list(factors=NULL, loadings=NULL),
				 nb_theta = 1, X = NULL, Z = NULL, sz = NULL){
  #Y is data with features=rows, observations=cols
  #L is number of desired latent dimensions
  #fam the likelihood for the data 
  #(poisson, negative binomial, binomial approximation to multinomial, bernoulli)
  #ctl a list of control parameters for optimization
  #penalty the L2 penalty for the latent factors
  #regression coefficients are not penalized
  #nb_theta see MASS::negative.binomial (nb_theta->infty = Poisson)
  #X a matrix of column (observations) covariates
  #any column with all same values (eg 1 for intercept) will be removed
  #this is because we force the intercept, so want to avoid collinearity
  #Z a matrix of row (feature) covariates, usually not needed
  #the basic model is R=AX'+ZG'+VU', where E[Y]=M=linkinv(R)
  #regression coefficients are A,G, latent factors are U and loadings V.
  
  fam<-match.arg(fam)
  N<-ncol(Y); J<-nrow(Y)
  #sanity check inputs
  if(fam %in% c("poi","nb","mult","bern")){ stopifnot(min(Y)>=0) }
  if(fam=="bern"){ stopifnot(max(Y)<=1) }
  
  #preprocess covariates and set updateable indices
  if(!is.null(X)){ 
    stopifnot(nrow(X)==ncol(Y))
    #we force an intercept, so remove it from X to prevent collinearity
    X<-remove_intercept(X)
    Ko<-ncol(X)+1
  } else {
    Ko<-1
  }
  if(!is.null(Z)){ 
    stopifnot(nrow(Z)==nrow(Y)) 
    Kf<-ncol(Z)
  } else {
    Kf<-0
  }
  lid<-(Ko+Kf)+(1:L)
  uid<-Ko + 1:(Kf+L)
  vid<-c(1:Ko, lid)
  Ku<-length(uid); Kv<-length(vid)
  
  #create glmpca_family object
  gnt<-glmpca_init(Y,fam,sz,nb_theta)
  gf<-gnt$gf; rfunc<-gnt$rfunc; a1<-gnt$intercepts
  
  #initialize U,V, with row-specific intercept terms
  U<-cbind(1, X, matrix(rnorm(N*Ku)*1e-5/Ku,nrow=N))
  if(!is.null(init$factors)){
    #print("initialize factors")
    L0<-min(L,ncol(init$factors))
    U[,(Ko+Kf)+(1:L0)]<-init$factors[,1:L0,drop=FALSE]
  }
  #a1 = naive MLE for gene intercept only
  V<-cbind(a1, matrix(rnorm(J*(Ko-1))*1e-5/Kv,nrow=J))
  V<-cbind(V, Z, matrix(rnorm(J*L)*1e-5/Kv,nrow=J))
  if(!is.null(init$loadings)){
    #print("initialize loadings")
    L0<-min(L,ncol(init$loadings))
    V[,(Ko+Kf)+(1:L0)]<-init$loadings[,1:L0,drop=FALSE]
  }
  
  #run optimization
  dev<-rep(NA,ctl$maxIter)
  for(t in 1:ctl$maxIter){
    #rmse[t]<-sd(Y-ilfunc(rfunc(U,V)))
    dev[t]<-gf$dev_func(Y,rfunc(U,V))
    if(t>5 && abs(dev[t]-dev[t-1])/(0.1+abs(dev[t-1]))<ctl$eps){
      break
    }
    if(verbose){ 
      dev_format<-format(dev[t],scientific=TRUE,digits=4)
      print(paste0("Iteration: ",t," | deviance=",dev_format)) 
    }
    #(k %in% lid) ensures no penalty on regression coefficients: 
    for(k in vid){
      ig<- gf$infograd(Y,rfunc(U,V))
      grads<- (ig$grad)%*%U[,k] - penalty*V[,k]*(k %in% lid) 
      infos<- (ig$info) %*% U[,k]^2 + penalty*(k %in% lid)
      V[,k]<-V[,k]+grads/infos
    }
    for(k in uid){
      ig<- gf$infograd(Y,rfunc(U,V))
      grads<- crossprod(ig$grad, V[,k]) - penalty*U[,k]*(k %in% lid) 
      infos<- crossprod(ig$info, V[,k]^2) + penalty*(k %in% lid) 
      U[,k]<-U[,k]+grads/infos
    }
  }
  #postprocessing: include row and columnn labels for regression coefficients
  if(is.null(Z)){
    G<-NULL
  } else {
    G<-U[,Ko+(1:Kf),drop=FALSE]
    rownames(G)<-colnames(Y); colnames(G)<-colnames(Z)
  }
  X<-if(is.null(X)){ matrix(1,nrow=N) } else { cbind(1,X) }
  if(!is.null(colnames(X))){ colnames(X)[1]<-"(Intercept)" }
  A<-V[,1:Ko,drop=FALSE]
  rownames(A)<-rownames(Y); colnames(A)<-colnames(X)
  res<-ortho(U[,lid],V[,lid],A,X=X,G=G,Z=Z,ret="df")
  rownames(res$factors)<-colnames(Y)
  rownames(res$loadings)<-rownames(Y)
  res$dev=dev[1:t]; res$family<-gf
  res
}
