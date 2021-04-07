# Warm Start Spatial Tweedie
pathMM_sptw_ridge <- function(y = NULL,
                              X = NULL,
                              Z = NULL,
                              index = NULL,
                              index.y.0 = NULL,
                              beta.init = NULL,
                              gamma.init = NULL,
                              alpha.init = NULL, 
                              lapMat=NULL,
                              p = 1.5,
                              tol = 1e-3,
                              miter = 1e4,
                              l1_seq=NULL,
                              l2_seq=NULL,
                              verbose=T){
  source('spatial_tweedie.R')
  # Initiate matrix for storage
  esteffb <- esteffa <- esteffg <- matrix(data=list(), nrow = length(l1_seq), ncol=1)
  nrows <- length(l1_seq)
  count <- 0
  for(i in 1:nrows){
    if(i==1){
      count <- count + 1
      print(c(i,count))
      pen_mat <- l1_seq[i]*diag(nrow(lapMat))
      
      temp_est <- spatial_tweedie(y = y,
                                  X = X,
                                  Z = Z,
                                  index = index,
                                  index.y.0 = index.y.0,
                                  beta.init = beta.init,
                                  gamma.init = gamma.init,
                                  alpha.init = alpha.init,
                                  pen_mat=pen_mat,
                                  p = p,
                                  tol = tol,
                                  miter = miter,
                                  verbose=verbose)
      esteffb[[i,1]] <- temp_est$beta
      esteffa[[i,1]] <- temp_est$alpha
      esteffg[[i,1]] <- temp_est$gamma
    }else{
      count <- count + 1
      print(c(i,count))
      betai <- esteffb[[(i-1),1]]
      alphai <- esteffa[[(i-1),1]]
      gammai <- esteffg[[(i-1),1]]
      
      pen_mat <- l1_seq[i]*diag(nrow(lapMat))
      temp_est <- spatial_tweedie(y = y,
                                  X = X,
                                  Z = Z,
                                  index = index,
                                  index.y.0 = index.y.0,
                                  beta.init = betai,
                                  gamma.init = gammai,
                                  alpha.init = alphai,
                                  pen_mat=pen_mat,
                                  p = p,
                                  tol = tol,
                                  miter = miter,
                                  verbose=verbose)
      esteffb[[i,1]] <- temp_est$beta
      esteffa[[i,1]] <- temp_est$alpha
      esteffg[[i,1]] <- temp_est$gamma
    }
  }
  list(beta=esteffb, alpha=esteffa, gamma=esteffg)
}
