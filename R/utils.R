#used by optimizers.R
stop_custom<-function(.subclass, message, call = NULL, ...){
  #https://adv-r.hadley.nz/conditions.html#signalling
  err <- structure(list(message=message,call=call,...),
                   class = c(.subclass, "error", "condition")
  )
  stop(err)
}

#used by optimizers.R
create_minibatches<-function(N,batch_size=1000,randomize=TRUE){
  #N = number of total observations in the data
  #if randomize=TRUE, we assume the data has not been shuffled,
  #...so minibatch indices are randomized
  #if randomize=FALSE, we assume the data itself 
  #...has been pre-shuffled so minibatches can be sequential
  batch_size<-max(batch_size,1) #batch size must be at least one observation
  batch_size<-min(batch_size,N) #batch size cannot be larger than dataset size
  B<-ceiling(N/batch_size) #round up to ensure no batch exceeds the requested batch_size
  i<-seq_len(N)
  if(B==1){ #case where only one minibatch containing whole dataset
    return(list(i))
  } else { #non-trivial number of minibatches
    bk<-cut(i,B,labels=FALSE)
    if(randomize) { #Y not preshuffled, must randomize minibatch indices
      return(split(sample.int(N,N,replace=FALSE),bk))
    } else { #Y is preshuffled, no need to randomize batches
      return(split(i,bk))
    }
  } 
}

#used by postprocess.R
l2norm<-function(x){sqrt(sum(x^2))}

#used by postprocess.R
colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  # apply(x,2,norm)
  sqrt(colSums(x^2))
}

# #used by initialization.R
# wlra<-function(Y,L,wts_row,wts_col,X=NULL,Z=NULL){
#   #represents the matrix Y with low rank approximation WH'
#   #where ncol(W)=ncol(H)=L and L is smaller than the dimensions of Y
#   #wts_row,wts_col are vectors of row and column weights
#   wr<-sqrt(wts_row); wc<-sqrt(wts_col)
#   Y<-wc*t(wr*Y) #this has dims of transposed Y
#   if(!is.null(X)){
#     fit<-lm.fit(X*wc,Y)
#     A<-t(coef(fit))/wr
#     Y<-t(residuals(fit))
#   } else {
#     A<-NULL
#     Y<-t(Y)
#   }
#   if(!is.null(Z)){
#     fit<-lm.fit(Z*wr,Y)
#     G<-t(coef(fit))/wc
#     Y<-residuals(fit)
#   } else {
#     G<-NULL
#   }
#   fit<-La.svd(Y, nu=L, nv=L)
#   V<-fit$u/wr
#   U<-t(fit$d[1:L]*fit$vt)/wc
#   list(A=A,G=G,U=U,V=V)
# }