#' @title GLM-PCA
#' @description Generalized principal components analysis for 
#'   dimension reduction of non-normally distributed data.
#' @name glmpca
#' 
#' @param Y matrix-like object of count or binary data with features as rows 
#'   and observations as columns. Sparse matrices from the \code{Matrix} 
#'   package are supported. Column-oriented sparsity is preferred.
#' @param L desired number of latent dimensions (positive integer).
#' @param fam string describing the likelihood to use for the data. Families
#'   include Poisson ('\code{poi}'), negative binomial with global 
#'   overdispersion ('\code{nb}'), negative binomial with feature-specific 
#'   overdispersion ('\code{nb2}'), or binomial ('\code{binom}'). Families 
#'   '\code{mult}' and '\code{bern}' are deprecated as both are special cases of
#'   '\code{binom}' with \code{sz} set to NULL and 1, respectively. They are 
#'   provided only for backward compatibility. Family '\code{nb2}' has not been
#'   thoroughly tested and is considered experimental.
#' @param minibatch string describing whether gradients should be computed with
#'   all observations ('\code{none}', the default) or a subset of observations, 
#'   which is useful for larger datasets. Option '\code{stochastic}' computes
#'   a noisy estimate of the full gradient using a random sample of observations
#'   at each iteration. Option '\code{memoized}' computes the full data 
#'   gradient under memory constraints by caching summary statistics across
#'   batches of observations.
#' @param optimizer string describing whether to use the fast AvaGrad method
#'   ('\code{avagrad}', the default) or the slower diagonal Fisher scoring 
#'   method ('\code{fisher}') that was used in the original glmpca 
#'   implementation.
#' @param ctl a list of control parameters. See 'Details'
#' @param sz numeric vector of size factors for each observation. If NULL 
#'   (the default), colSums are used for family '\code{binom}', and 
#'   colMeans are used for families '\code{poi}','\code{nb}', and '\code{nb2}'.
#' @param nb_theta initial value for negative binomial overdispersion 
#'   parameter(s). Small values lead to more overdispersion. Default: 100. See
#'   \code{\link[MASS]{negative.binomial}}. (\code{nb_theta}->\eqn{\infty}
#'   equivalent to Poisson).
#' @param X a matrix of column (observations) covariates. Any column with all
#'   same values (eg. 1 for intercept) will be removed. This is because we force
#'   a feature-specific intercept and want to avoid collinearity.
#' @param Z a matrix of row (feature) covariates, usually not needed.
#' @param init a list containing initial estimates for the factors (\code{U}) 
#'   and loadings (\code{V}) matrices.
#' @param ... additional named arguments. Provided only for backward 
#'   compatibility.
#'
#' @details The basic model is \eqn{R = AX'+ZG'+VU'}, where \eqn{E[Y] = M
#'   = linkinv(R)}. Regression coefficients are \code{A} and \code{G}, latent
#'   factors are \code{U} and loadings are \code{V}. 
#'   The objective is to minimize the deviance between \code{Y} 
#'   and \code{M}. The deviance quantifies the goodness-of-fit of the GLM-PCA
#'   model to the data (smaller=better). 
#'   Note that \code{glmpca} uses a random initialization, 
#'   so for fully reproducible results one may use \code{set.seed}.
#'   
#'   The \code{ctl} argument accepts any of the following optional components:
#'   \describe{
#'     \item{verbose}{Logical. Should detailed status messages be printed 
#'       during the optimization run? Default: \code{FALSE}.}
#'     \item{batch_size}{Positive integer. How many observations should be
#'       included in a minibatch? Larger values use more memory but lead to 
#'       more accurate gradient estimation. Ignored if \code{minibatch='none'}.
#'       Default: 1000.}
#'     \item{lr}{Positive scalar. The AvaGrad learning rate. Large values
#'       enable faster convergence but can lead to numerical instability.
#'       Default: 0.1. If a numerical divergence occurs, \code{glmpca}
#'       will restart the optimization \code{maxTry} times (see below)
#'       and reduce the learning rate by a factor of five each time.}
#'     \item{penalty}{Non-negative scalar. The L2 penalty for the latent 
#'       factors. Default: 1. Regression coefficients are not penalized. Only
#'       used by the Fisher scoring optimizer. Larger values improve numerical
#'       stability but bias the parameter estimates. If a numerical divergence 
#'       occurs, \code{glmpca} will restart the optimization \code{maxTry} times
#'       (see below) and increase the penalty by a factor of five each time.}
#'     \item{maxTry}{Positive integer. In case of numerical divergence, how
#'       many times should optimization be restarted with a more stable penalty
#'       or learning rate? Default: 10.}
#'     \item{minIter}{Positive integer. Minimum number of iterations (full
#'       passes through the dataset) before checking for numerical convergence.
#'       Default: 30.}
#'     \item{maxIter}{Positive integer. Maximum number of iterations. If
#'       numerical convergence is not achieved by this point, the results may
#'       not be reliable and a warning is issued. Default: 1000.}
#'     \item{tol}{Positive scalar. Relative tolerance for assessing convergence.
#'       Convergence is determined by comparing the deviance at the previous
#'       iteration to the current iteration. Default: 1e-4.}
#'     \item{epsilon}{Positive scalar. Avagrad hyperparameter. See Savarese et 
#'       al (2020). Default: 0.1.}
#'     \item{betas}{Numeric vector of length two. Avagrad hyperparameters. 
#'       See Savarese et al (2020). Default: \code{c(0.9, 0.999)}.}
#'     \item{minDev}{Scalar. Minimum deviance threshold at which optimization
#'       is terminated. Useful for comparing different algorithms as it avoids
#'       the need to determine numerical convergence. Default: NULL}
#'   }
#'   
#' @return An S3 object of class \code{glmpca} with copies of input components 
#'   \code{optimizer}, \code{minibatch}, \code{ctl},\code{X}, and \code{Z},
#'   along with the following additional fitted components:
#'   \describe{
#'     \item{factors}{a matrix \code{U} whose rows match the columns 
#'       (observations) of \code{Y}. It is analogous to the principal components
#'       in PCA. Each column of the factors matrix is a different latent 
#'       dimension.}
#'     \item{loadings}{a matrix \code{V} whose rows match the rows 
#'       (features/dimensions) of \code{Y}. It is analogous to loadings in PCA. 
#'       Each column of the loadings matrix is a different latent dimension.}
#'     \item{coefX}{a matrix \code{A} of coefficients for the 
#'       observation-specific covariates matrix \code{X}. Each row of coefX 
#'       corresponds to a row of \code{Y} and each column corresponds to a 
#'       column of \code{X}. The first column of coefX contains feature-specific 
#'       intercepts which are included by default.}
#'     \item{coefZ}{a matrix \code{G} of coefficients for the feature-specific 
#'       covariates matrix \code{Z}. Each row of coefZ corresponds to a column 
#'       of \code{Y} and each column corresponds to a column of \code{Z}. By 
#'       default no such covariates are included and this is returned as NULL.}
#'     \item{dev}{a vector of deviance values. The length of the vector is the 
#'       number of iterations it took for GLM-PCA's optimizer to converge. 
#'       The deviance should generally decrease over time. 
#'       If it fluctuates wildly, this often indicates numerical instability, 
#'       which can be improved by decreasing the learning rate or increasing the 
#'       penalty, see \code{ctl}.}
#'     \item{dev_smooth}{a locally smoothed version of \code{dev} that may be
#'       easier to visualize when \code{minibatch='stochastic'}.}
#'     \item{glmpca_family}{an S3 object of class glmpca_family. This is a minor
#'       extension to the \link[stats]{family} or \link[MASS]{negative.binomial}
#'       object used by functions like \link[stats]{glm} and 
#'       \link[MASS]{glm.nb}. It is basically a list with various internal 
#'       functions and parameters needed to optimize the GLM-PCA objective 
#'       function. For the negative binomial case, it also contains the final 
#'       estimated value of the overdispersion parameter (\code{nb_theta}).}
#'     \item{offsets}{For Poisson and negative binomial families, the offsets
#'       are the logarithmically transformed size factors. These are needed to
#'       compute the predicted mean values.}
#'   }
#' 
#' @examples
#' #create a simple dataset with two clusters
#' mu<-rep(c(.5,3),each=10)
#' mu<-matrix(exp(rnorm(100*20)),nrow=100)
#' mu[,1:10]<-mu[,1:10]*exp(rnorm(100))
#' clust<-rep(c("red","black"),each=10)
#' Y<-matrix(rpois(prod(dim(mu)),mu),nrow=nrow(mu))
#' #visualize the latent structure
#' res<-glmpca(Y, 2)
#' factors<-res$factors
#' plot(factors[,1],factors[,2],col=clust,pch=19)
#' 
#' @seealso
#' \code{\link{predict.glmpca}}, 
#' \code{\link[stats]{prcomp}}, \code{\link[stats]{glm}},
#' \code{\link[logisticPCA]{logisticSVD}},
#' \code{\link[scry]{devianceFeatureSelection}}, 
#' \code{\link[scry]{nullResiduals}}
#' 
#' @references
#' Savarese P, McAllester D, Babu S, and Maire M (2020).
#' Domain-independent Dominance of Adaptive Methods. \emph{arXiv}
#' \url{https://arxiv.org/abs/1912.01823}
#' 
#' Townes FW (2019). Generalized Principal Component Analysis. \emph{arXiv}
#' \url{https://arxiv.org/abs/1907.02647}
#' 
#' Townes FW, Hicks SC, Aryee MJ, and Irizarry RA (2019).
#' Feature Selection and Dimension Reduction for Single Cell RNA-Seq based on a 
#' Multinomial Model. \emph{Genome Biology}
#' \url{https://doi.org/10.1186/s13059-019-1861-6}
#' 
#' @importFrom methods is
#' @export
glmpca<-function(Y, L, fam=c("poi","nb","nb2","binom","mult","bern"), 
                   minibatch=c("none","stochastic","memoized"),
                   optimizer=c("avagrad","fisher"), ctl = list(), 
                   sz=NULL, nb_theta=NULL, X=NULL, Z=NULL, 
                   init=list(factors=NULL, loadings=NULL), ...){
  #Y is a matrix-like object, must support the following operations:
  #min,max,as.matrix,sum,colSums,colMeans,rowSums,rowMeans,`[`
  
  #this line is purely for backward compatibility with scry v1.0.0 (bioc3.11)
  if(length(fam)>1){ fam<-fam[1] }
  fam<-match.arg(fam)
  minibatch<-match.arg(minibatch)
  optimizer<-match.arg(optimizer)
  #handle deprecated arguments from old function signature
  if(fam=="mult"){
    message("Family 'mult' is deprecated. Coercing to equivalent family ",
            "'binom' with sz=NULL instead. Please use this in the future, ",
            "as 'mult' will eventually be removed.")
    fam<-"binom"; sz<-NULL
  }
  if(fam=="bern"){
    message("Family 'bern' is deprecated. Coercing to equivalent family ",
            "'binom' with sz=1 instead. Please use this in the future, ",
            "as 'bern' will eventually be removed.")
    fam<-"binom"; sz<-1
  }
  dots<-list(...)
  if(!is.null(dots$penalty)){
    message("Control parameter 'penalty' should be provided as an element ",
            "of 'ctl' rather than a separate argument.")
    if(is.null(ctl$penalty)){ ctl$penalty<-dots$penalty }
  }
  if(!is.null(dots$verbose)){
    message("Control parameter 'verbose' should be provided as an element ",
            "of 'ctl' rather than a separate argument.")
    if(is.null(ctl$verbose)){ ctl$verbose<-dots$verbose }
  }
  
  N<-ncol(Y); J<-nrow(Y)
  #sanity check inputs
  if(fam %in% c("poi","nb","nb2","binom")){ stopifnot(min(Y) >= 0) }
  if(!is.null(sz)){
    stopifnot(all(sz>0))
    if(fam=="binom"){
      #really this should be all(colMax(Y)<=sz), will fix later since rare
      stopifnot(max(Y) <= max(sz)) 
    }
  }
  
  ic<-init_ctl(N,fam,minibatch,optimizer,ctl)
  fam<-ic$fam; minibatch<-ic$minibatch; optimizer<-ic$optimizer; ctl<-ic$ctl
  if(minibatch=="none"){ 
    if(is(Y,"sparseMatrix")){
      message("Sparse matrices are not supported for minibatch='none'. ",
              "Coercing to dense matrix. If this exhausts memory, ",
              "consider setting minibatch to 'stochastic' or 'memoized'.")
    }
    Y<-as.matrix(Y)
  } #in future, handle DelayedArray too
  
  #create glmpca_family object
  gnt<-glmpca_init(Y, fam, sz=sz, nb_theta=nb_theta)
  gf<-gnt$gf; rfunc<-gnt$rfunc; offsets<-gnt$offsets
  
  #initialize factors and loadings matrices
  uv<-uv_init(N, J, L, gnt$intercepts, X=X, Z=Z, init=init)
  U<-uv$U; V<-uv$V; lid<-uv$lid; uid<-uv$uid; vid<-uv$vid
  
  #minimize the deviance using an optimizer
  fit<-NULL
  for(ntry in seq.int(ctl$maxTry)){
    if(optimizer=="avagrad"){ 
      if(ctl$verbose){
        message("Trying AvaGrad with learning rate: ",signif(ctl$lr,4))
      }
      e<-tryCatch(
        if(minibatch=="none"){
          fit<-avagrad_optimizer(Y,U,V,uid,vid,ctl,gf,rfunc,offsets)
        } else if(minibatch=="stochastic"){
          fit<-avagrad_stochastic_optimizer(Y,U,V,uid,vid,ctl,gf,rfunc,offsets)
        } else {
          fit<-avagrad_memoized_optimizer(Y,U,V,uid,vid,ctl,gf,rfunc,offsets)
        },
        error_glmpca_divergence=function(e){ e },
        error_glmpca_dev_incr=function(e){ e }
      )
      if(is.null(fit)){ ctl$lr<-ctl$lr/5 } else { break }
    } else if(optimizer=="fisher"){
      if(ctl$verbose){
        message("Trying Fisher scoring with penalty: ",ctl$penalty)
      }
      e<-tryCatch(
        fit<-fisher_optimizer(Y,U,V,uid,vid,ctl,gf,rfunc,offsets),
        error_glmpca_divergence=function(e){ e },
        error_glmpca_dev_incr=function(e){ e }
      )
      if(is.null(fit)){ 
        if(ctl$penalty==0){ 
          ctl$penalty<-1
        } else { 
          ctl$penalty<-ctl$penalty*5 
        }
      } else { 
        break 
      }
      # } else if(optimizer=="none"){ #closed-form approximate likelihood only
      #   dev<-gf$dev_func(Y,rfunc(U,V,offsets))
      #   fit<-list(U=U, V=V, dev=dev, gf=gf)
      #   break
    } else {
      stop("unrecognized optimizer")
    }
  } #end for loop over different penalty or learning rate values
  if(ntry==ctl$maxTry){ stop(e) }
  if(length(fit$dev)==ctl$maxIter){
    warning("Reached maximum number of iterations (",ctl$maxIter,
            ") without numerical convergence. Results may be unreliable.")
  }
  #print(class(Y))
  res<-postprocess(fit,uid,vid,lid,rnames=rownames(Y),cnames=colnames(Y))
  fpars<-list(ctl=ctl,offsets=offsets,optimizer=optimizer,minibatch=minibatch)
  res<-c(res,fpars) #S3 object of class "glmpca" (really just a big list)
  class(res)<-"glmpca"
  res
}

#' @importFrom utils tail
#' @export
print.glmpca<-function(x,...){
  cat("GLM-PCA fit with", ncol(x$factors), "latent factors")
  cat("\nnumber of observations:",nrow(x$factors))
  cat("\nnumber of features:",nrow(x$loadings))
  cat("\nfamily:",x$glmpca_family$glmpca_fam)
  cat("\nminibatch:",x$minibatch)
  cat("\noptimizer:",x$optimizer)
  if(x$optimizer=="fisher"){
    cat("\nl2 penalty:",x$ctl$penalty)
  }
  if(x$optimizer=="avagrad"){
    cat("\nlearning rate:",signif(x$ctl$lr,3))
  }
  dev_final<-format(tail(x$dev,1),scientific=TRUE,digits=4)
  cat("\nfinal deviance:",dev_final)
  invisible(x)
}

#' @title Predict Method for GLM-PCA Fits
#' @description Predict the mean matrix from a fitted generalized principal
#'   component analysis model object.
#' @name predict.glmpca
#' 
#' @param object a fitted object of class inheriting from \code{glmpca}.
#' @param ... additional named arguments. Currently ignored.
#' 
#' @details Let \code{Y} be the data matrix originally used to estimate the
#'   parameters in \code{fit}. The GLM-PCA model regards each element of 
#'   \code{Y} as a random sample from an exponential family distribution 
#'   such as a Poisson, negative binomial, or binomial likelihood. The 
#'   components of a GLM-PCA fit are combined to produce the predicted 
#'   mean of this distribution at each entry of \code{Y}. This matrix may be
#'   regarded as a 'denoised' version of the original data.
#' 
#' @return a dense \code{matrix} of predicted mean values.
#' 
#' @section Warning:
#'   The predicted mean matrix returned by this function
#'   will have the same dimensions as the original data matrix and it will be
#'   dense even if the original data were sparse. This can exhaust available
#'   memory for large datasets, so use with caution.
#'   
#' @seealso 
#' \code{\link{glmpca}}, 
#' \code{\link[stats]{predict.glm}} with \code{type='response'}
#' 
#' @export
predict.glmpca<-function(object,...){#output="matrix"){
  #given a fitted glmpca object, return the fitted mean matrix
  #can be useful for imputation of noisy data
  #warning: will produce a dense matrix which can overwhelm memory for large
  #datasets
  gf<-object$glmpca_family
  offsets<-object$offsets
  binom_n<-gf$binom_n
  ilfunc<-gf$linkinv
  if(is.null(object$coefZ)){
    U<-as.matrix(cbind(object$X, object$factors))
    V<-as.matrix(cbind(object$coefX, object$loadings))
  } else {
    U<-as.matrix(cbind(object$X, object$coefZ, object$factors))
    V<-as.matrix(cbind(object$coefX, object$Z, object$loadings))
  }
  if(is.null(offsets) || all(offsets==0)){ #everything but poi, nb, nb2
    M<-ilfunc(tcrossprod(V,U))
  } else { #poi, nb, nb2
    M<-ilfunc(t(offsets+tcrossprod(U,V)))
  }
  if(is.null(binom_n) || all(binom_n)==1){ #everything but binomial w/ n>=2
    return(M)
  } else if(length(binom_n)>1) { #binomial approx to multinomial
    return(t(t(M)*binom_n))
  } else { #binomial w/ global n>=2 (eg n=2 for SNP data
    return(M*binom_n)
  }
}
