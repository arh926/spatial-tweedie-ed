crossvalPll_sptw <- function(K=3,
                             y=NULL,
                             X=NULL,
                             Z=NULL,
                             index=NULL,
                             beta.init = NULL,
                             gamma.init = NULL,
                             alpha.init = NULL,
                             id.list=NULL,
                             l1_seq=NULL,
                             l2_seq=NULL,
                             lapMat=NULL,
                             miter=NULL,
                             tol=NULL,
                             p=NULL,
                             verbose=T){
  nrows <- length(l1_seq); ncols <- length(l2_seq)
  
  require(foreach)
  require(doParallel)
  registerDoParallel()
  
  nCores = detectCores()
  
  if(K<=nCores){
    cl = makeCluster(K)
    registerDoParallel(cl, cores = K)
  } else{
    cl = makeCluster(nCores)
    registerDoParallel(cl, cores = nCores)
  }
  
  crossval.dev = foreach(i=1:K) %dopar% {
    crossvalMat <- matrix(NA, nrow=nrows, ncol=ncols)
    require(tweedie)
    source('pathMM_sptw.R')
    source("crossvalPll_sptw.R")
    # print sink file
    if(verbose) sink(paste("mdfold_exact",i,".txt", sep=""))
    
    id.train <- as.numeric(unlist(id.list[-i]))
    estEff <- pathMM_sptw(y = y[id.train],
                          X = X[id.train,],
                          Z = Z[id.train,],
                          index = index[id.train],
                          index.y.0 = (y[id.train]==0),
                          beta.init = beta.init,
                          gamma.init = gamma.init,
                          alpha.init = alpha.init, 
                          lapMat=lapMat,
                          p = p,
                          tol = tol,
                          miter = miter,
                          l1_seq=l1_seq,
                          l2_seq=l2_seq,
                          verbose=verbose)
    
    id.test <- as.numeric(unlist(id.list[i]))
    effectsi.beta <- estEff$beta
    effectsi.alpha <- estEff$alpha
    effectsi.gamma <- estEff$gamma
    estp <- estEff$p
    y.test <- y[id.test]
    X.test <- X[id.test,]
    index.test <- index[id.test]
    
    # unlink(paste("mdfold_exact",i,".txt", sep=""))
    for(k in 1:ncols){
      for(l in 1:nrows){
        beta.eff <- effectsi.beta[l,k][[1]]
        alpha.eff <- effectsi.alpha[l,k][[1]]
        p <- estp[l,k]
        mu.test <- exp(crossprod(t(X.test),beta.eff)+alpha.eff[match(index.test,names(alpha.eff))])
        
        crossvalMat[l,k] <- sum(tweedie.dev(y=y.test,mu=mu.test,power = p))
      }
    }
    list(dev=crossvalMat,eff.beta=effectsi.beta, eff.alpha=effectsi.alpha, eff.gamma=effectsi.gamma, eff.p=estp)
  }
  on.exit(stopCluster(cl));gc(TRUE)
  crossval.dev
}
