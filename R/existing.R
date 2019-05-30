# Functions for running existing methods such as tSNE, PCA, ZIFA, etc on censored data
library(Matrix)

rm_zero_rowcol<-function(Y){
  #remove all rows and columns containing all zeros
  Y<-Y[rowSums(Y>0)>0,] #remove rows with zeros all the way across
  Y<-Y[,colSums(Y>0)>0]
  Y
}

norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}

pca_irlba<-function(Y,L=2,center=TRUE,scale=TRUE,rmzero=TRUE){
  #Y is data matrix, L is desired number of latent dimensions
  #scale is flag of whether or not to center and scale Y before applying irlba
  #Y<-as.matrix(Y)
  if(rmzero==TRUE) Y<-rm_zero_rowcol(Y)
  #if(scale) Y<-t(scale(t(Y)))
  #svd1<-irlba::irlba(Y,L)
  #factors<-t(svd1$d * t(svd1$v)) #using recycling
  factors<-irlba::prcomp_irlba(t(Y),n=L,center=center,scale=scale)
  #colnames(factors)<-paste0("pca_irlba",1:L)
  colnames(factors)<-paste0("dim",1:L)
  as.data.frame(factors$x)
}

pca<-function(Y,L=2,center=TRUE,scale=TRUE,rmzero=TRUE,ret_obj=FALSE){
  Y<-as.matrix(Y)
  if(rmzero==TRUE) Y<-rm_zero_rowcol(Y)
  res<-prcomp(as.matrix(t(Y)),center=center,scale.=scale,rank.=L)
  factors<-as.data.frame(res$x)
  colnames(factors)<-paste0("dim",1:L)
  if(ret_obj){
    return(list(factors=factors,obj=res))
  } else{
    return(factors)
  }
}

mds<-function(Y,L=2,metric=TRUE,distance="euclidean",scale=TRUE,rmzero=TRUE){
  #Multidimensional Scaling
  #Y is data
  #L is desired latent dimension
  #metric=TRUE means use cmdscale(). metric=FALSE means use isoMDS()
  #see http://www.statmethods.net/advstats/mds.html
  #distance is choice of distance function passed to dist()
  if(rmzero==TRUE) Y<-rm_zero_rowcol(Y)
  Yt<-scale(t(Y),scale=scale)
  d <- dist(as.matrix(Yt),distance) # euclidean distances between the cols
  if(metric){
    fit<-cmdscale(d,k=L)
  } else {
    fit<-MASS::isoMDS(d,k=L)$points
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

tsne<-function(Y,L=2,center=TRUE,scale=TRUE,rmzero=TRUE,method="Rtsne",...){
  if(rmzero==TRUE){
    Yt<-t(rm_zero_rowcol(Y))
  } else {
    Yt<-t(Y)
  }
  if(center || scale){
    Yt<-scale(Yt,center=center,scale=scale)
  }
  Yt<-as.matrix(Yt)
  if(method=="Rtsne"){
    fit<-Rtsne::Rtsne(Yt,dims=L,...)$Y
  } else if(method=="tsne"){
    fit<-tsne::tsne(Yt,k=L,...)
  }
  colnames(fit)<-paste0("dim",1:L)
  as.data.frame(fit)
}

simlr<-function(Y,L=2,...){
  Y<-as.matrix(Y)
  nclust<-SIMLR::SIMLR_Estimate_Number_of_Clusters(Y)
  nclust<-which.min(nclust$K1+nclust$K2)
  #suppress printing of verbose SIMLR messages
  simlr_output<-capture.output(res<-SIMLR::SIMLR(Y,nclust,no.dim=L,...))
  #the clustering results are in res$y
  factors<-as.data.frame(res$ydata)
  colnames(factors)<-paste0("dim",1:L)
  factors
}

cidr<-function(Y,L=2){
  #assumes Y is unnormalized counts not log transformed
  sData <- cidr::scDataConstructor(as.matrix(Y))
  sData <- cidr::determineDropoutCandidates(sData)
  sData <- cidr::wThreshold(sData)
  sData <- cidr::scDissim(sData)
  sData <- cidr::scPCA(sData, plotPC=FALSE)
  sData <- cidr::nPC(sData)
  factors<-as.data.frame(sData@PC[,1:L])
  colnames(factors)<-paste0("dim",1:L)
  factors
}

zinbwave<-function(Y,L=2,parallel=FALSE){
  #Y is unnormalized counts not log transformed
  if(parallel){
    bp<-BiocParallel::bpparam()
  } else {
    bp<-BiocParallel::SerialParam()
  }
  suppressWarnings(fit<-zinbwave::zinbFit(as.matrix(Y), K=L, BPPARAM=bp))
  factors<-as.data.frame(zinbwave::getW(fit))
  colnames(factors)<-paste0("dim",1:L)
  rownames(factors)<-colnames(Y)
  factors
}

ppca<-function(Y,L=2,...){
  stop("implementation not yet finished")
  pcaMethods::pp(Y,L,...)
}

regpca<-function(Y,L=2,X=NULL,center=TRUE,scale=TRUE,rmzero=TRUE){
  #regress against covariates X
  #then do PCA on residuals
  #if X is NULL, will use detection rates and intercept as covars
  #... is for additional arguments passed to function 'pca'
  #... for example center, scale, rmzero
  Y<-as.matrix(Y)
  if(rmzero==TRUE) Y<-rm_zero_rowcol(Y)
  if(is.null(X)) X<-matrix(colMeans(Y>0)) #detection rates for each cell
  detlm<-lm(t(Y)~X) #multivariate linear regression
  Yrs<-t(resid(detlm))
  #Ysvd<-svd(Yrs)
  ##V<-t(Ysvd$u)[1:L,]
  #factors<-as.data.frame(t(Ysvd$d[1:L]*t(Ysvd$v)[1:L,]))
  #colnames(factors)<-paste0("dim",1:L)
  #factors
  pca(Yrs,L=L,center=center,scale=scale,rmzero=FALSE)
}