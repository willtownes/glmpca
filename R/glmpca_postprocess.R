# rotate GLM-PCA factors to orthonormal basis

l2norm<-function(x){sqrt(sum(x^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  # apply(x,2,norm)
  sqrt(colSums(x^2))
}

ortho<-function(U, V, A, X=rep(1,nrow(U)), ret=c("m","df"), G=NULL, Z=NULL, 
                rnames=NULL, cnames=NULL){
  #U is NxL matrix of cell factors, V is JxL matrix of loadings onto genes
  #X is NxKo matrix of cell specific covariates
  #A is JxKo matrix of coefficients of X
  #Z is JxKf matrix of gene specific covariates
  #G is NxKf matrix of coefficients of Z
  #ret= return either data frame or matrix format
  #imputed expression: E[Y] = g^{-1}(R) where R = AX'+ZG'+VU'
  #rnames=row names of original data matrix
  #cnames=column names of original data matrix
  
  ret<-match.arg(ret)
  L<-ncol(U)
  X<-as.matrix(X,nrow=nrow(U))
  #if(is.null(Z)){ Z<-0 }
  if(is.null(Z) || all(Z==0) || all(G==0)){
    G<-Z<-NULL
  }
  if(!is.null(G)){
    if(length(Z)==1){ Z<-rep(Z,nrow(V)) }
    Z<-as.matrix(Z,nrow=nrow(V))
  }
  #we assume A is not null or zero
  #remove correlation between U and A
  #at minimum, this will cause factors to have mean zero
  reg<-lm.fit(X,U)
  factors<-residuals(reg)
  A<-A+tcrossprod(V,coef(reg))
  colnames(A)<-colnames(X)
  #remove correlation between V and G
  if(is.null(G)){
    loadings<-V
  } else { #G is not empty
    reg<-lm.fit(Z,V)
    loadings<-residuals(reg)
    G<-G+tcrossprod(factors,coef(reg))
    colnames(G)<-colnames(Z)
  }
  #rotate factors to make loadings orthornormal
  svdres<-svd(loadings)
  loadings<-svdres$u
  factors<-t(t(factors%*%svdres$v)*svdres$d)
  #arrange latent dimensions in decreasing L2 norm
  o<-order(colNorms(factors),decreasing=TRUE)
  factors<-factors[,o,drop=FALSE]
  loadings<-loadings[,o,drop=FALSE]
  colnames(loadings)<-colnames(factors)<-paste0("dim",1:L)
  if(!is.null(cnames)){
    rownames(factors)<-rownames(X)<-cnames
    if(!is.null(G)){ rownames(G)<-cnames }
  }
  if(!is.null(rnames)){
    rownames(loadings)<-rownames(A)<-rnames
    if(!is.null(Z)){ rownames(Z)<-rnames }
  }
  if(ret=="df"){
    loadings<-as.data.frame(loadings)
    factors<-as.data.frame(factors)
    X<-as.data.frame(X)
    A<-as.data.frame(A)
    if(!is.null(G)){ 
      Z<-as.data.frame(Z)
      G<-as.data.frame(G)
    }
  }
  list(factors=factors, loadings=loadings, X=X, Z=Z, coefX=A, coefZ=G)
}

postprocess<-function(fit,uid,vid,lid,rnames=NULL,cnames=NULL){
  #fit is a list containing optimized params U,V, the dev trace, and gf object
  #postprocessing: include row and column labels for regression coefficients
  X<-fit$U[,-uid,drop=FALSE] #will always have at least one column for intercepts
  if(!is.null(colnames(X))){ colnames(X)[1]<-"(Intercept)" }
  A<-fit$V[,-uid,drop=FALSE]
  if(length(vid)==ncol(fit$V)){ 
    Z<-G<-NULL
  } else if(length(vid)<ncol(fit$V)) { 
    Z<-fit$V[,-vid,drop=FALSE] 
    G<-fit$U[,-vid,drop=FALSE]
  } else { stop("length(vid) cannot exceed ncol(fit$V)") }
  res<-ortho(fit$U[,lid,drop=FALSE], fit$V[,lid,drop=FALSE], A, X=X, 
             ret="df", G=G, Z=Z, rnames=rnames, cnames=cnames)
  res$dev<-fit$dev; res$glmpca_family<-fit$gf
  res$dev_smooth<-fit$dev_smooth #NULL unless minibatch=='stochastic'
  res
}
