# Warm Start Spatial Tweedie
pathMM_sptw <- function(y = NULL,
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
  esteffb <- esteffa <- esteffg <- matrix(data=list(), nrow = length(l1_seq), ncol = length(l2_seq))
  p_est <- matrix(NA, nrow = length(l1_seq), ncol = length(l2_seq))
  nrows <- length(l1_seq); ncols <- length(l2_seq)
  count <- 0
  for(i in 1:nrows){
    if(i%%2==1){
      for(j in 1:ncols){
        count <- count + 1
        print(c(i,j,count))
        if(j==1 & i==1){
          pen_mat <- l1_seq[i]*diag(nrow(lapMat))+l2_seq[j]*lapMat
          
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
          esteffb[[i,j]] <- temp_est$beta
          esteffa[[i,j]] <- temp_est$alpha
          esteffg[[i,j]] <- temp_est$gamma
          p_est[i,j] <- temp_est$p
        }else if(j==1 & i!=1){
          betai <- esteffb[[(i-1),1]]
          alphai <- esteffa[[(i-1),1]]
          gammai <- esteffg[[(i-1),1]]
          pi <- p_est[(i-1),1]
          
          pen_mat <- l1_seq[i]*diag(nrow(lapMat))+l2_seq[j]*lapMat
          temp_est <- spatial_tweedie(y = y,
                                      X = X,
                                      Z = Z,
                                      index = index,
                                      index.y.0 = index.y.0,
                                      beta.init = betai,
                                      gamma.init = gammai,
                                      alpha.init = alphai,
                                      pen_mat=pen_mat,
                                      p = pi,
                                      tol = tol,
                                      miter = miter,
                                      verbose=verbose)
          esteffb[[i,j]] <- temp_est$beta
          esteffa[[i,j]] <- temp_est$alpha
          esteffg[[i,j]] <- temp_est$gamma
          p_est[i,j] <- temp_est$p
        }else{
          betai <- esteffb[[i,(j-1)]]
          alphai <- esteffa[[i,(j-1)]]
          gammai <- esteffg[[i,(j-1)]]
          pi <- p_est[i,(j-1)]
          
          pen_mat <- l1_seq[i]*diag(nrow(lapMat))+l2_seq[j]*lapMat
          temp_est <- spatial_tweedie(y = y,
                                      X = X,
                                      Z = Z,
                                      index = index,
                                      index.y.0 = index.y.0,
                                      beta.init = betai,
                                      gamma.init = gammai,
                                      alpha.init = alphai,
                                      pen_mat=pen_mat,
                                      p = pi,
                                      tol = tol,
                                      miter = miter,
                                      verbose=verbose)
          esteffb[[i,j]] <- temp_est$beta
          esteffa[[i,j]] <- temp_est$alpha
          esteffg[[i,j]] <- temp_est$gamma
          p_est[i,j] <- temp_est$p
        }
      }
    }else{
      for(j in ncols:1){
        count <- count + 1
        print(c(i,j,count))
        if(j==ncols){
          betai <- esteffb[[(i-1),j]]
          alphai <- esteffa[[(i-1),j]]
          gammai <- esteffg[[(i-1),j]]
          pi <- p_est[(i-1),j]
          
          pen_mat <- l1_seq[i]*diag(nrow(lapMat))+l2_seq[j]*lapMat
          temp_est <- spatial_tweedie(y = y,
                                      X = X,
                                      Z = Z,
                                      index = index,
                                      index.y.0 = index.y.0,
                                      beta.init = betai,
                                      gamma.init = gammai,
                                      alpha.init = alphai,
                                      pen_mat=pen_mat,
                                      p = pi,
                                      tol = tol,
                                      miter = miter,
                                      verbose=verbose)
          esteffb[[i,j]] <- temp_est$beta
          esteffa[[i,j]] <- temp_est$alpha
          esteffg[[i,j]] <- temp_est$gamma
          p_est[i,j] <- temp_est$p
        }else{
          betai <- esteffb[[i,(j+1)]]
          alphai <- esteffa[[i,(j+1)]]
          gammai <- esteffg[[i,(j+1)]]
          pi <- p_est[i,(j+1)]
          
          pen_mat <- l1_seq[i]*diag(nrow(lapMat))+l2_seq[j]*lapMat
          temp_est <- spatial_tweedie(y = y,
                                      X = X,
                                      Z = Z,
                                      index = index,
                                      index.y.0 = index.y.0,
                                      beta.init = betai,
                                      gamma.init = gammai,
                                      alpha.init = alphai,
                                      pen_mat=pen_mat,
                                      p = pi,
                                      tol = tol,
                                      miter = miter,
                                      verbose=verbose)
          esteffb[[i,j]] <- temp_est$beta
          esteffa[[i,j]] <- temp_est$alpha
          esteffg[[i,j]] <- temp_est$gamma
          p_est[i,j] <- temp_est$p
        }
      }
    }
  }
  list(beta=esteffb, alpha=esteffa, gamma=esteffg, p=p_est)
}
