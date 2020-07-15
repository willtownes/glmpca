#An object of class glmpca_family is very similar to the standard 'family' objects used by glm()
#The differences are:
# 1. glmpca_family does not have the 'aic', 'initialize', or 'simfun' components that are in family objects.
# 2. glmpca_family has additional components 'glmpca_fam', 'infograd', and 'dev_func'
# ...glmpca_fam is a string identifying the likelihood for the data (eg 'poi', 'nb' etc)
# ...infograd is a function of data matrix Y and real-valued linear predictor R (often called eta in GLM literature)
# ...infograd returns a list with 'grad' (the gradient with respect to R) and 'info' (the Fisher information)
# ...note infograd is for maximizing the log-likelihood. To minimize deviance, take the negative of both terms.
# ...dev_func is also a function of data Y and real-valued predictor R, it returns the total deviance between the 
# ...fitted model implied by R and the glmpca_family (containing any dispersion parameters) and a saturated model where mean=data.
# 3. glmpca_family supports negative binomial likelihood where theta can be a vector instead of just a scalar.
# 4. For poisson and negative binomial families, the offsets vector is stored under 'offsets'
# 4. For negative binomial family, the theta vector is stored under 'nb_theta'
# 5. For binomial family, the size factors vector is stored under 'binom_n'. It can also be a scalar.

print.glmpca_family <- function(x, ...){
  #based on print.family eg 
  #https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/family.R#L21
  cat("\nGLM-PCA fam:", x$glmpca_fam, "\n")
  cat("Family:", x$family, "\n")
  cat("Link function:", x$link, "\n\n")
  invisible(x)
}

#' @importFrom stats make.link
nb2_family<-function(theta=stop("'theta' must be specified")){
  #modification of MASS::negative.binomial family that allows theta to be a vector instead of a scalar
  #only log link function is supported.
  #we assume this will be applied to data with observations in columns and features in rows
  #ie the vector theta aligns with the rows of the data matrix Y and the mean matrix Mu.
  #this assumption facilitates recycling in the variance function etc.
  linktemp <- "log"
  stats<-make.link(linktemp)
  ilfunc <- stats$linkinv
  variance <- function(Mu){ Mu + Mu^2/theta } #recycling
  validmu <- function(Mu){ all(Mu > 0) }
  dev.resids<-function(Y, Mu, Wt=1){ #note: if wt is provided it must be a matrix of same dims as Y and Mu
    2*Wt*( Y*log(pmax(1,Y)/Mu) - (Y+theta)*log((Y+theta)/(Mu+theta)))
  }
  infograd<-function(Y,R,sz=NULL,grad_only=TRUE){
    M<-ilfunc(R) #ilfunc=exp
    W<-1/variance(M)
    list(grad=(Y-M)*W*M, info=if(grad_only){ NULL } else { W*M^2 })
  }
  dev_func<-function(Y,R,sz=NULL){
    #Note dev.resids function gives the square of actual residuals
    #so the sum of these is the deviance
    sum(dev.resids(Y,ilfunc(R),1)) 
  }
  famname <- paste0("Negative binomial (theta vector length ",length(theta),")")
  #stuff that is also found in 'family' objects from base R
  gf<-list(family = famname, link = linktemp, linkfun = stats$linkfun, 
           linkinv = ilfunc, variance = variance, 
           dev.resids = dev.resids, mu.eta = stats$mu.eta, 
           validmu = validmu, valideta = stats$valideta)
  #add the GLM-PCA specific stuff
  gf2<-list(glmpca_fam="nb2", nb_theta=theta, 
            infograd=infograd, dev_func=dev_func)
  structure(c(gf,gf2), class="glmpca_family")
}

mat_binom_dev<-function(X,P,n){
  #binomial deviance for two matrices
  #X,P are JxN matrices
  #n is vector of length N (same as cols of X,P) or a scalar
  #nz<-X>0
  
  #ensure recycling works properly
  if(length(n)>1){ X<-t(X); P<-t(P) }
  term1<-sum(X*log(X/(n*P)), na.rm=TRUE)
  #nn<-x<n
  nx<-n-X
  term2<-sum(nx*log(nx/(n*(1-P))), na.rm=TRUE)
  2*(term1+term2)
}

#' @importFrom stats binomial poisson
glmpca_family<-function(fam,binom_n=NULL,nb_theta=NULL){
  #binom_n a scalar or vector of total count params for binomial distribution
  #if binom_n is a vector, it is a binomial approx to multinomial (n=total count for each observation)
  #global binom_n=1 is useful for binary (0/1) data (ie Bernoulli distr)
  #global binom_n=2 may be useful for SNP data.
  #nb_theta a scalar or vector of overdispersion params (smaller theta=more overdispersion)
  #if nb_theta is a vector it is one value per feature
  #if nb_theta is a scalar it is a global overdispersion shared by all features
  
  #create GLM family object
  if(fam %in% c("poi","nb","nb2")){
    if(fam=="poi"){
      gf<-poisson(link="log")
    } else if(fam=="nb"){
      gf<-MASS::negative.binomial(nb_theta, link="log")
      gf$nb_theta<-nb_theta
    } else if(fam=="nb2"){
      return(nb2_family(nb_theta)) #everything already handled by nb2_family function
    }
    #gf$offsets<-gf$linkfun(sz)
  } else if(fam=="binom"){
    gf<-binomial(link="logit")
    if(is.null(binom_n)){ stop("size factor(s) 'binom_n' must be provided if fam='binom'") }
    gf$binom_n<-binom_n
  } else {
    stop("unrecognized family type")
  }
  #remove family components we don't need for GLM-PCA
  gf$aic <- gf$initialize <- gf$simfun <- NULL
  #change class and add glmpca_fam component
  class(gf)<-"glmpca_family"
  gf$glmpca_fam<-fam
  #variance function, determined by GLM family
  vfunc<-gf$variance
  #inverse link func, mu as a function of linear predictor R
  ilfunc<-gf$linkinv
  if(fam=="poi"){
    gf$infograd<-function(Y,R,sz=NULL,grad_only=TRUE){
      M<-ilfunc(R) #ilfunc=exp
      list(grad=(Y-M), info=if(grad_only){ NULL } else { M })
    }
  } else if(fam=="nb"){
    gf$infograd<-function(Y,R,sz=NULL,grad_only=TRUE){
      M<-ilfunc(R) #ilfunc=exp
      W<-1/vfunc(M)
      list(grad=(Y-M)*W*M, info=if(grad_only){ NULL } else { W*M^2 })
    }
  } else if(fam=="binom"){
    if(length(binom_n)>1){
      gf$infograd<-function(Y,R,sz,grad_only=TRUE){
        Pt<-t(ilfunc(R)) #ilfunc=expit, Pt may be very small probabilities
        list(grad=Y-t(sz*Pt), 
             info=if(grad_only){ NULL } else { t(sz*vfunc(Pt))})
      }
    } else if(binom_n==1){ #Bernoulli
      gf$infograd<-function(Y,R,sz=1,grad_only=TRUE){
        P<-ilfunc(R)
        list(grad=Y-P, info=if(grad_only){ NULL} else { vfunc(P) })
      }
    } else { #binomial with global N, eg N=2 for modeling SNPs 
      gf$infograd<-function(Y,R,sz=binom_n,grad_only=TRUE){
        P<-ilfunc(R)
        list(grad=Y-sz*P, info=if(grad_only){ NULL} else { sz*vfunc(P) })
      }
    }
  } else { #this is not actually used but keeping for future reference
    #this is most generic formula for GLM but computationally slow
    stop("invalid fam")
    #derivative of inverse link function, dmu/dR
    hfunc<-gf$mu.eta
    gf$infograd<-function(Y,R,sz=NULL,grad_only=TRUE){
      M<-ilfunc(R)
      W<-1/vfunc(M)
      H<-hfunc(R) #hfunc(R)
      list(grad=(Y-M)*W*H, info=if(grad_only){ NULL } else { W*H^2 })
    }
  }
  #create deviance function
  if(fam=="binom" && any(binom_n!=1)){
    gf$dev_func<-function(Y,R,sz){
      mat_binom_dev(Y,ilfunc(R),sz)
    }
  } else {
    gf$dev_func<-function(Y,R,sz=NULL){
      #Note dev.resids function gives the square of actual residuals
      #so the sum of these is the deviance
      sum(gf$dev.resids(Y,ilfunc(R),1)) 
    }
  }
  gf
}
