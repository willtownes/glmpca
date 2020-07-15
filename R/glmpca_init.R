#initialization functions for GLM-PCA
#source("./algs/glmpca2/glmpca_families.R")
#source("./algs/glmpca2/function_approximation.R")

init_ctl<-function(N,fam,minibatch,optimizer,ctl){
  #set default control params
  if(!is.null(ctl$eps)){
    message("Control parameter 'eps' is deprecated. Coercing to equivalent ",
            "'tol'. Please use 'tol' in the future as 'eps' will eventually ",
            "be removed")
    ctl$tol<-ctl$eps
  }
  if(is.null(ctl$verbose)){ ctl$verbose<-FALSE }
  if(is.null(ctl$minIter)){ 
    ctl$minIter<-20 
  } else {
    stopifnot(ctl$minIter>=1)
  }
  if(is.null(ctl$maxIter)){ 
    ctl$maxIter<-1000 
  } else { 
    stopifnot(ctl$maxIter>ctl$minIter)
  }
  if(is.null(ctl$maxTry)){ 
    ctl$maxTry<-10 
  } else { 
    stopifnot(ctl$maxTry>=1) 
  }
  if(is.null(ctl$tol)){ 
    #ctl$tol<-if(optimizer=="memoized"){ 1e-5 } else { 1e-4 }
    ctl$tol<-1e-4
  } else {
    stopifnot(ctl$tol>0)
  }
  #minibatch defaults
  if(minibatch != "none"){
    if(is.null(ctl$batch_size)){ 
      ctl$batch_size<-1000 
    } else {
      stopifnot(ctl$batch_size>=1)
    }
    if(ctl$batch_size>=N){ 
      message("batch_size exceeds number of data points, ",
              "setting minibatch to 'none'.")
      minibatch<-"none" 
    }
  }
  if(minibatch=="none"){ #Y is dense matrix, no minibatches
    #Y<-as.matrix(Y)
    ctl$batch_size<-NULL
  } else { #either stochastic gradient or memoization
    if(optimizer!="avagrad"){ 
      msg<-"forcing optimizer to 'avagrad' to enable minibatch="
      message(msg,minibatch) 
      optimizer<-"avagrad"
    }
    if(fam=="nb2"){
      msg1<-"Family 'nb2' is only supported for minibatch=none, "
      msg2<-"using 'nb' instead."
      message(msg1,msg2)
      fam<-"nb"
    }
    # if(minibatch=="memoized" && fam=="nb"){
    #   msg<-"Family 'nb' is not yet supported for memoized minibatches."
    #   stop(msg)
    # }
  } 
  #optimizer defaults
  if(optimizer=="fisher"){
    if(is.null(ctl$penalty)){ 
      ctl$penalty<-1 
    } else {
      stopifnot(ctl$penalty>=0)
    }
  } else {
    ctl$penalty<-NULL
  }
  if(optimizer=="avagrad"){
    if(is.null(ctl$lr)){ ctl$lr<-0.1 } else { stopifnot(ctl$lr>0) }
    if(is.null(ctl$epsilon)){ 
      ctl$epsilon<-0.1 
    } else { 
      stopifnot(ctl$epsilon>0) 
    }
    if(is.null(ctl$betas)){ 
      ctl$betas<-c(.9,.999) 
    } else {
      stopifnot(length(ctl$betas)==2)
      stopifnot(min(ctl$betas)>=0)
      stopifnot(max(ctl$betas)<1)
    }
  } else {
    ctl$lr<-ctl$epsilon<-ctl$betas<-NULL
  }
  #minDev is also a control parameter initialized to NULL by default
  list(fam=fam,minibatch=minibatch,optimizer=optimizer,ctl=ctl)
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
  #if sz is NULL, the default depends on 'fam'
  #poi,nb,nb2: colMeans(Y)
  #binom: colSums(Y) (approximates a multinomial)
  #note: binomial with sz=1 (ie Bernoulli) may be useful for binary data
  #note: binomial with sz=2 may be useful for SNP data
  N<-ncol(Y); J<-nrow(Y)
  #check that Matrix package is available
  if(is(Y,"Matrix")){
    # http://r-pkgs.had.co.nz/description.html#dependencies
    if(!requireNamespace("Matrix",quietly=TRUE)){
      stop("package 'Matrix' required when Y is a Matrix. Please install it.")
    } else {
      colSums<-Matrix::colSums; rowSums<-Matrix::rowSums
      colMeans<-Matrix::colMeans; rowMeans<-Matrix::rowMeans
    }
  } # else if(is(Y,"DelayedArray")){
    # if(!requireNamespace("DelayedArray",quietly=TRUE)){
    #   stop("package 'DelayedArray' required when Y is a DelayedArray. Please install it.")
    # } else {
    #   colSums<-DelayedArray::colSums; rowSums<-DelayedArray::rowSums
    #   colMeans<-DelayedArray::colMeans; rowMeans<-DelayedArray::rowMeans
    # }
  # }
  has_offset<-fam %in% c("poi","nb","nb2")
  #sz must be either a scalar or one value per column of Y
  if(is.null(sz)){ 
    if(has_offset){
      sz<-colMeans(Y)
    } else if(fam=="binom"){
      sz<-colSums(Y)
    } else {
      sz<-1
    }
  } else {
    stopifnot(length(sz) %in% c(1,ncol(Y))) 
  }
  binom_n<-offsets<-NULL
  if(has_offset){
    offsets<-log(sz)
  } else if(fam=="binom"){
    binom_n<-sz
  } #else fam=="gaussian" both stay NULL
  if(is.null(nb_theta)){
    nb_theta<-switch(fam, nb=100.0, nb2=rep(100.0,J))
  }
  gf<-glmpca_family(fam,binom_n=binom_n,nb_theta=nb_theta)
  
  if(all(sz==1)){ #Bernoulli and Gaussian
    a1<-gf$linkfun(rowMeans(Y))
  } else if(length(sz)>1) { # Poisson, nb, binomial approx to multinom
    a1<-gf$linkfun(rowSums(Y)/sum(sz))
  } else { #single global size factor, eg Binomial(2,p) for SNP data
    a1<-gf$linkfun(rowSums(Y)/(N*sz))
  }
  if(any(is.infinite(a1))){
    stop("Some rows were all zero, please remove them.")
  }
  #rfunc is the real-valued linear predictor
  if(has_offset){
    rfunc<-function(U,V,offsets){ t(offsets+tcrossprod(U,V)) } 
  } else { #no offsets, eg fam 'binom' or gaussian
    #size factors are incorporated via binom fam object
    rfunc<-function(U,V,offsets=NULL){ tcrossprod(V,U) }
  }
  list(gf=gf, rfunc=rfunc, intercepts=a1, offsets=offsets)
}

wlra<-function(Y,L,wts_row,wts_col,X=NULL,Z=NULL){
  #represents the matrix Y with low rank approximation WH'
  #where ncol(W)=ncol(H)=L and L is smaller than the dimensions of Y
  #wts_row,wts_col are vectors of row and column weights
  wr<-sqrt(wts_row); wc<-sqrt(wts_col)
  Y<-wc*t(wr*Y) #this has dims of transposed Y
  if(!is.null(X)){
    fit<-lm.fit(X*wc,Y)
    A<-t(coef(fit))/wr
    Y<-t(residuals(fit))
  } else {
    A<-NULL
    Y<-t(Y)
  }
  if(!is.null(Z)){
    fit<-lm.fit(Z*wr,Y)
    G<-t(coef(fit))/wc
    Y<-residuals(fit)
  } else {
    G<-NULL
  }
  fit<-La.svd(Y, nu=L, nv=L)
  V<-fit$u/wr
  U<-t(fit$d[1:L]*fit$vt)/wc
  list(A=A,G=G,U=U,V=V)
}

uv_init<-function(N, J, L, a1, X=NULL, Z=NULL, 
                  init=list(factors=NULL, loadings=NULL)){
  #N: number of observations (columns) of data matrix
  #J: number of features (rows) of data matrix
  #L=number of desired latent dimensions
  #a1=feature-specific intercept term (something like the log(rowMeans) of Y)
  #X=covariates matched to observations of Y (cols)
  #Z=covariates matched to features of Y (rows)
  #init: a list with following elements (each is optional)
  #factors: initialized values of observation-specific factors
  #loadings: initialized values of feature-specific loadings

  #preprocess covariates and set updateable indices
  if(is.null(X)){ 
    X<-matrix(1,nrow=N,ncol=1)
  } else {
    stopifnot(nrow(X)==N)
    #we force an intercept, so remove it from X to prevent collinearity
    X<-remove_intercept(X)
    X<-cbind(1,X)
  }
  Ko<-ncol(X)
  if(is.null(Z)){ 
    Kf<-0
  } else {
    stopifnot(nrow(Z)==J) 
    Kf<-ncol(Z)
  }
  lid<-(Ko+Kf)+(1:L)
  uid<-Ko + 1:(Kf+L)
  vid<-c(1:Ko, lid)
  #Ku<-Kf+L; Kv<-Ko+L
  Ku<-length(uid); Kv<-length(vid)
  U<-cbind(X, matrix(rnorm(N*Ku,sd=1e-5/Ku),nrow=N))
  #initialize U,V, with row-specific intercept terms
  if(!is.null(init$factors)){
    L0<-min(L,ncol(init$factors))
    U[,(Ko+Kf)+(1:L0)]<-init$factors[,1:L0,drop=FALSE]
  }
  #a1 = naive MLE for gene intercept only
  V<-cbind(a1, matrix(rnorm(J*(Ko-1),sd=1e-5/Kv),nrow=J)) #=A matrix
  V<-cbind(V, Z, matrix(rnorm(J*L,sd=1e-5/Kv),nrow=J)) #=A,Z,V matrices
  if(!is.null(init$loadings)){
    L0<-min(L,ncol(init$loadings))
    V[,(Ko+Kf)+(1:L0)]<-init$loadings[,1:L0,drop=FALSE]
  }
  # } else if(init=="pal"){
  #   if(!is.matrix(Y) && N>1000){
  #     ss<-sample.int(N,1000)
  #     Y<-Y[,ss]
  #     if(!is.null(offsets)){ offsets<-offsets[ss] }
  #     if(!is.null(gf$binom_n)){ gf$binom_n<-gf$binom_n[ss] }
  #   }
  #   Y<-as.matrix(Y)
  #   fam<-gf$glmpca_fam
  #   if(fam %in% c("poi","nb","nb2")){
  #     #quadratic approximation to Poisson likelihood 
  #     #followed by closed-form weighted low-rank approximation (WLRA)
  #     sz<-exp(offsets)
  #     f<-function(l,h){quad_approx(exp,l,h)}
  #     cba<-t(mapply(f,a1-4,a1+4)) #matrix with cols c,b,a
  #   } else if(fam=="binom"){
  #     #quadratic approximation to binomial likelihood 
  #     #followed by closed-form weighted low-rank approximation (WLRA)
  #     sz<-gf$binom_n
  #     f<-function(l,h){quad_approx(softplus,l,h)}
  #     cba<-t(mapply(f,a1-6,a1+6))
  #   } else {
  #     stop("unrecognized likelihood")
  #   }
  #   #this is the part that creates a dense matrix in-memory
  #   Y<-(Y-outer(cba[,"b"],sz))/(2*outer(cba[,"a"],sz))
  #   fit<-wlra(Y,L,cba[,"a"],sz,X=X,Z=Z)
  #   #U<-cbind(X, fit$G, fit$U)
  #   #U is already initialized randomly, discard wlra estimates since ss of cells
  #   V<-cbind(fit$A, Z, fit$V)
  # } else { 
  #   stop("unrecognized 'init' value")
  # }
  list(U=U,V=V,lid=lid,uid=uid,vid=vid)
}
