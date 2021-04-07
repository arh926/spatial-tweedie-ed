require(tweedie)
require(Matrix)
# Simulation
########################################################
# generate data: covariates + spatial effect: response #
########################################################
# setwd("/Users/aritrahalder/Dropbox/PhD Research/paper 1/code/paper1_code/") # setwd("/Users/aritrahalder/")
# setwd("P:/Desktop/paper1_code")
# setwd("/Users/aritrahalder/Desktop/paper1_code/")

###############
# Spatial     #
###############
# require(sp)
# require(spdep)
# require(rgeos)
# require(rgdal)
# require(sf)
# require(MBA)
# generate spatial effect
# blockwise
# pdf("/Users/aritrahalder/Desktop/speff.pdf")


# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],speff),ct_shape)
# plot(ct_shape)
# points(ct_zip[,c("longitude","latitude")], pch="+", col=col,cex=.5)


# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],speff),ct_shape)
# plot(ct_shape)
# points(ct_zip[,c("longitude","latitude")], pch="+", col=col,cex=.5)


# dev.off()

###############
# Non-Spatial #
###############

# mu <- y; phi <- phi_sim; p <- 1.5; epsilon <- 1e-16
## Not Run:: Just Check if we have the same likelihood
# index.y.0 <- y==0
# series_likelihood <- function(y=NULL, p=NULL, mu=NULL, phi=NULL, epsilon=1e-16, index.y.0=NULL){
#   alpha <- (2-p)/(p-1)
#   -sum(1/phi*(y*mu^(1-p)/(1-p)-mu^(2-p)/(2-p)))-sum(log(apply(cbind(y[!index.y.0],phi[!index.y.0]),1, function(x){
#       k=x[1]^(2-p)/((2-p)*x[2])
#       w_max <- (x[1]/(p-1))^(k*alpha)/(2-p)^k/factorial(k)/factorial(k*alpha-1)/x[2]^(k*(1+alpha))
#       k_min <- 1
#       while((x[1]/(p-1))^(k_min*alpha)/(2-p)^k_min/factorial(k_min)/factorial(k_min*alpha-1)/x[2]^(k_min*(1+alpha))<epsilon*w_max) k_min <- k_min+1
#       k_max <- k+1
#       while((x[1]/(p-1))^(k_max*alpha)/(2-p)^k_max/factorial(k_max)/factorial(k_max*alpha-1)/x[2]^(k_max*(1+alpha))>epsilon*w_max) k_max <- k_max+1
#       return(sum((x[1]/(p-1))^((k_min:k_max)*alpha)/(2-p)^(k_min:k_max)/factorial((k_min:k_max))/factorial((k_min:k_max)*alpha-1)/x[2]^((k_min:k_max)*(1+alpha)))/x[1])
#     })))
# }
# Test:: series_likelihood(y,1.5, mu_sim,phi_sim, index.y.0=index.y.0)
# -sum(log(dtweedie.series(y,1.5, mu_sim,phi_sim)))


# saddle_likelihood <- function(y=NULL, p=NULL, mu=NULL, phi=NULL, epsilon=1/6){
#  lik_vec <- 1/phi*((y^(2-p)-y*mu^(1-p))/(1-p)-(y^(2-p)-mu^(2-p))/(2-p))+log(2*pi*phi*(y+epsilon)^p)/2
#  lik_vec[y==0] <- mu[y==0]^(2-p)/(phi[y==0]*(2-p))
#  sum(lik_vec)
# }
# -sum(log(dtweedie.saddle(y,1.5, mu_sim,phi_sim)))
# saddle_likelihood(y,1.5, mu_sim,phi_sim)




####################
# Cross Validation #
####################
# split into testing and training
N <- c(1e4,2e4,3e4,5e4)
psiL=c(7,10,27,110)
psiU=c(10,27,110,240)

setwd("/Users/aritrah/OneDrive - University of Connecticut/PhD Research/paper 1/latex+code/code")
load("adjCT.RData")
adjM = Matrix(adjCT)
degM = Matrix(diag(rowSums(adjCT))); rownames(degM) <- colnames(degM) <- rownames(adjCT)
source('spatial_tweedie.R'); source('sp_plot-1.R'); source('crossvalPll_sptw.R'); source('fold_split.R')
source('crossvalPll_sptw_ridge.R'); source('pathMM_sptw_ridge.R')
nrep <- 10
beta <-  beta_mle <- beta_ridge <- array(NA,c(8,nrep,4, 4))
alpha <- alpha_mle <- alpha_ridge <- array(NA,c(nrow(adjCT), nrep,4, 4))
gamma <- gamma_mle <- gamma_ridge <- array(NA,c(8,nrep,4, 4))
lambda_est <- array(NA,c(2,nrep,4, 4)); lambda_est_ridge <- array(NA,c(2,nrep,4, 4))
for(i in 1:length(N)){
  for(j in 1:length(psiL)){
    for(k in 1:nrep){
      load(paste("/Users/aritrah/OneDrive - University of Connecticut/PhD Research/paper 1/latex+code/code/sim_data/blockwise_",i,j,".RData",sep=""))
      sim_data <- obj$data
      x_mat <- obj$X
      z_mat <- obj$Z
      p <- 0.6
      zipnames <- rownames(adjCT)
      train <- unlist(sapply(zipnames, function(x){
        id <- which(sim_data$zip==x)
        train_id <- sample(id, round(p*length(id)), replace = F)
        train_id
      }))
      train_test_ind <- rep("InT_V",N[i])
      train_test_ind[train] <- "InT_T"
      sim_data$train_ind <- train_test_ind
      data_train <- sim_data[sim_data$train_ind=="InT_T",]
      
      full_id <- fold_split(K=3,index = data_train$zip)
      fold_1 <- as.numeric(unlist(lapply(full_id, function(x) x[[1]])))
      fold_2 <- as.numeric(unlist(lapply(full_id, function(x) x[[2]])))
      fold_3 <- as.numeric(unlist(lapply(full_id, function(x) x[[3]])))
      #fold_4 <- as.numeric(unlist(lapply(full_id, function(x) x[[4]])))
      #fold_5 <- as.numeric(unlist(lapply(full_id, function(x) x[[5]])))
      
      id.list <- list(fold1=fold_1, fold2=fold_2, fold3=fold_3)#,fold4=fold_4,fold5=fold_5)
      
      y.train = data_train$y
      X.train = x_mat[sim_data$train_ind=="InT_T",]
      Z.train = z_mat[sim_data$train_ind=="InT_T",]
      index = data_train$zip
      beta.init = rep(0,8)#runif(8)
      gamma.init = rep(0,8)#runif(8)
      alpha.init = rep(0,nrow(adjCT)); names(alpha.init) <- rownames(adjCT)
     
      p = 1.2
      tol = 1e-6
      miter = 1e4
      l1_seq <- exp(seq(-5,5,length.out = 10))
      l2_seq <- exp(seq(-5,5,length.out = 10))
      lapMat <- degM - adjM
      
      cvM <- crossvalPll_sptw(K=3,
                              y=y.train,
                              X=X.train,
                              Z=Z.train,
                              index=index,
                              beta.init = beta.init,
                              gamma.init = gamma.init,
                              alpha.init = alpha.init,
                              id.list=id.list,
                              l1_seq=l1_seq,
                              l2_seq=l2_seq,
                              lapMat=lapMat,
                              miter=miter,
                              tol=tol,
                              p=p,
                              verbose=T)
      #installr::kill_all_Rscript_s()
      devM <- (cvM[[1]]$dev+cvM[[2]]$dev+cvM[[3]]$dev)
      pM <- (cvM[[1]]$eff.p+cvM[[2]]$eff.p+cvM[[3]]$eff.p)/3
      lambda_est[,k,j,i] <- arr.min <- which(devM==min(devM),arr.ind = T)
      
      # Ridge
      cvV <- crossvalPll_sptw_ridge(K=3,
                                    y=y.train,
                                    X=X.train,
                                    Z=Z.train,
                                    index=index,
                                    beta.init = beta.init,
                                    gamma.init = gamma.init,
                                    alpha.init = alpha.init,
                                    id.list=id.list,
                                    l1_seq=l1_seq,
                                    l2_seq=l2_seq,
                                    lapMat=lapMat,
                                    miter=miter,
                                    tol=tol,
                                    p=p,
                                    verbose=F)
      #installr::kill_all_Rscript_s()
      devV <- (cvV[[1]]$dev+cvV[[2]]$dev+cvV[[3]]$dev)
      lambda_est_ridge[,k,j,i] <- arr.minv <- which(devV==min(devV),arr.ind = T)
      
      # MLE
      fit_mle <- spatial_tweedie(y = y.train,
                                 X = X.train,
                                 Z = Z.train,
                                 index = index,
                                 index.y.0 = y.train==0,
                                 beta.init = beta.init,
                                 gamma.init = gamma.init,
                                 alpha.init = alpha.init,
                                 pen_mat = 0*diag(nrow(lapMat))+0*lapMat, #
                                 p = 1.8,
                                 tol = tol,
                                 miter = miter,
                                 inf=T)
      # Ridge
      fit_ridge <- spatial_tweedie(y = y.train,
                                   X = X.train,
                                   Z = Z.train,
                                   index = index,
                                   index.y.0 = y.train==0,
                                   beta.init = beta.init,
                                   gamma.init = gamma.init,
                                   alpha.init = alpha.init,
                                   pen_mat = l1_seq[arr.minv[1]]*diag(nrow(lapMat)),
                                   p = p,
                                   tol = tol,
                                   miter = miter,
                                   inf=T)
      
      # Proposed
      fit_sptw <- spatial_tweedie(y = y.train,
                                  X = X.train,
                                  Z = Z.train,
                                  index = index,
                                  index.y.0 = y.train==0,
                                  beta.init = beta.init,
                                  gamma.init = gamma.init,
                                  alpha.init = alpha.init,
                                  pen_mat = l1_seq[arr.min[1]]*diag(nrow(lapMat))+l2_seq[arr.min[2]]*lapMat, #
                                  p = pM[arr.min[1],arr.min[2]],
                                  tol = tol,
                                  miter = miter,
                                  inf=T,
                                  p.update = F)
      beta[,k,j,i] <- fit_sptw$optim_pars$beta
      alpha[,k,j,i] <- fit_sptw$optim_pars$alpha
      gamma[,k,j,i] <- fit_sptw$optim_pars$gamma
      
      beta_ridge[,k,j,i] <- fit_ridge$optim_pars$beta
      alpha_ridge[,k,j,i] <- fit_ridge$optim_pars$alpha
      gamma_ridge[,k,j,i] <- fit_ridge$optim_pars$gamma
      
      beta_mle[,k,j,i] <- fit_mle$optim_pars$beta
      alpha_mle[,k,j,i] <- fit_mle$optim_pars$alpha
      gamma_mle[,k,j,i] <- fit_mle$optim_pars$gamma
      
      print(c(i,j,k))
    }
  }
}

################
# Diagnostics  #
################
betaTrue <- obj$betaT; gammaTrue=obj$gammaT; speff=obj$alphaT
# pdf("diagnostic_blockwise.pdf")
plot3D::persp3D(x = log(l1_seq),
                y = log(l2_seq),
                z = pM,
                phi = 50,
                theta = 40,
                alpha=0.6,
                contour=T,
                expand=1.1,
                colkey=list(side=1),
                horizontal=T,
                bty="u",
                facets=T)
plot3D::persp3D(x = log(l1_seq),
                y = log(l2_seq),
                z = devM,
                phi = 50,
                theta = -40,
                alpha=0.6,
                contour=T,
                expand=1.1,
                colkey=list(side=1),
                horizontal=T,
                bty="u",
                facets=T)
plot(x=log(l1_seq), y=devV, type="l")

par(mfcol=c(1,3))
niter <- nrow(fit_mle$optim_pars$objFn)
obj_full <- rep(NA, 2*niter)
for(i in 1:niter){
  obj_full[2*i-1] <- fit_mle$optim_pars$objFn[i,1]
  obj_full[2*i] <- fit_mle$optim_pars$objFn[i,2]
} 
plot(obj_full, type="c", ylab="Objective Function", xlab="Iteration", lwd=2)
points(cbind(seq(1,2*niter, by=2),obj_full[seq(1,2*niter, by=2)]), col="darkgreen",  pch="*", cex=1.5)
points(cbind(seq(2,2*niter, by=2),obj_full[seq(2,2*niter, by=2)]), col="darkred",pch="o", cex=1.2)
legend("topright", c("Mean update", "p and Dispersion update"), pch=c("*","o"), col=c("darkgreen","darkred"))
grid()

niter <- nrow(fit_ridge$optim_pars$objFn)
obj_full <- rep(NA, 2*niter)
for(i in 1:niter){
  obj_full[2*i-1] <- fit_ridge$optim_pars$objFn[i,1]
  obj_full[2*i] <- fit_ridge$optim_pars$objFn[i,2]
} 
plot(obj_full, type="c", ylab="Objective Function", xlab="Iteration", lwd=2)
points(cbind(seq(1,2*niter, by=2),obj_full[seq(1,2*niter, by=2)]), col="darkgreen",  pch="*", cex=1.5)
points(cbind(seq(2,2*niter, by=2),obj_full[seq(2,2*niter, by=2)]), col="darkred", cex=1.5)
grid()

niter <- nrow(fit_sptw$optim_pars$objFn)
obj_full <- rep(NA, 2*niter)
for(i in 1:niter){
  obj_full[2*i-1] <- fit_sptw$optim_pars$objFn[i,1]
  obj_full[2*i] <- fit_sptw$optim_pars$objFn[i,2]
} 
plot(obj_full, type="c", ylab="Objective Function", xlab="Iteration", lwd=2)
points(cbind(seq(1,2*niter, by=2),obj_full[seq(1,2*niter, by=2)]), col="darkgreen",  pch="*", cex=1.5)
points(cbind(seq(2,2*niter, by=2),obj_full[seq(2,2*niter, by=2)]), col="darkred", cex=1.5)
grid()
par(mfcol=c(1,1))
cbind(mle=as.vector(c(round(sum((fit_mle$optim_pars$alpha-speff)^2),3),
                      round(cor(fit_mle$optim_pars$alpha,speff),3),
                      round(sum((fit_mle$optim_pars$beta-betaTrue)^2),3),
                      round(sum((fit_mle$optim_pars$gamma-gammaTrue)^2),3))),
      ridge=as.vector(c(round(sum((fit_ridge$optim_pars$alpha-speff)^2),3),
                        round(cor(fit_ridge$optim_pars$alpha,speff),3),
                        round(sum((fit_ridge$optim_pars$beta-betaTrue)^2),3),
                        round(sum((fit_ridge$optim_pars$gamma-gammaTrue)^2),3))),
      proposed=as.vector(c(round(sum((fit_sptw$optim_pars$alpha-speff)^2),3),
                  round(cor(fit_sptw$optim_pars$alpha,speff),3),
                  round(sum((fit_sptw$optim_pars$beta-betaTrue)^2),3),
                  round(sum((fit_sptw$optim_pars$gamma-gammaTrue)^2),3))))


# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],fit_sptw$optim_pars$alpha),ct_shape)
# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],speff),ct_shape)
par(mfcol=c(1,3))
plot(x=fit_mle$optim_pars$alpha,y=speff,
     ylab="True",
     xlab="Estimated",)
abline(c(0,1), col="darkred",lwd=2); grid()
plot(x=fit_ridge$optim_pars$alpha,y=speff,
     ylab="True",
     xlab="Estimated",)
abline(c(0,1), col="darkred",lwd=2); grid()
plot(x=fit_sptw$optim_pars$alpha,y=speff,
     ylab="True",
     xlab="Estimated",)
abline(c(0,1), col="darkred",lwd=2); grid()

boxplot(fit_mle$optim_pars$alpha~round(speff,2),
        ylab="True",
        xlab="Estimated",ylim=c(-5,3)); abline(h=unique(speff), col="darkred",lwd=2)
boxplot(fit_ridge$optim_pars$alpha~round(speff,2),
        ylab="True",
        xlab="Estimated",ylim=c(-5,3)); abline(h=unique(speff), col="darkred",lwd=2)
boxplot(fit_sptw$optim_pars$alpha~round(speff,2),
        ylab="True",
        xlab="Estimated",ylim=c(-5,3)); abline(h=unique(speff), col="darkred",lwd=2)
par(mfcol=c(1,1))
#plot(x=fit_sptw$optim_pars$beta,y=betaTrue, ylab="true (beta)",xlab="estimated (beta)");abline(c(0,1),col="red");grid()
#plot(x=fit_sptw$optim_pars$gamma,y=gammaTrue, ylab="true (gamma)",xlab="estimated (gamma)");abline(c(0,1),col="red");grid()
# dev.off()
##############
# Prediction #
##############
data_test <- sim_data[sim_data$train_ind=="InT_V",]
X.test = x_mat[sim_data$train_ind=="InT_V",]
Z.test = z_mat[sim_data$train_ind=="InT_V",]

# deviance ratio
beta_est <- fit_sptw$optim_pars$beta
alpha_est <- fit_sptw$optim_pars$alpha
Xbest <- crossprod(t(X.test),beta_est)+alpha_est[match(data_test$zip,names(alpha_est))]
Xbtrue <- crossprod(t(X.test),betaTrue)+data_test$speff


dev_est <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbest),power = 1.5))
dev_true <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbtrue),power = 1.5))
dev_est/dev_true

# deviance ratio
beta_est <- fit_ridge$optim_pars$beta
alpha_est <- fit_ridge$optim_pars$alpha
Xbest <- crossprod(t(X.test),beta_est)+alpha_est[match(data_test$zip,names(alpha_est))]
Xbtrue <- crossprod(t(X.test),betaTrue)+data_test$speff


dev_est <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbest),power = 1.5))
dev_true <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbtrue),power = 1.5))
dev_est/dev_true

# deviance ratio
beta_est <- fit_mle$optim_pars$beta
alpha_est <- fit_mle$optim_pars$alpha
Xbest <- crossprod(t(X.test),beta_est)+alpha_est[match(data_test$zip,names(alpha_est))]
Xbtrue <- crossprod(t(X.test),betaTrue)+data_test$speff


dev_est <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbest),power = 1.5))
dev_true <- sum(tweedie.dev(y = data_test$y,mu = exp(Xbtrue),power = 1.5))
dev_est/dev_true


