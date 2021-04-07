setwd("/Users/aritrahalder/Desktop/paper1_code/")
rm(list=ls())
require(tweedie)
require(zipcode)
data(zipcode)
ct_zip <- subset(zipcode,state=="CT")
# us.shp <- st_read("/Users/aritrahalder/Desktop/tl_2018_us_zcta510/tl_2018_us_zcta510.shp")
# ct.shp <- us.shp[!is.na(match(as.character(us.shp$ZCTA5CE10),ct_zip$zip)),]

usa_shape <- readRDS("GADM_2.8_USA_adm2.rds")
ct_shape <- subset(usa_shape,NAME_1=="Connecticut")
# plot(ct_shape)
# points(ct_zip[,c("longitude","latitude")], pch="+", col="blue",cex=.5)

# Not Run: Create Adjacency Matrix
# adjCT <- nb2mat(poly2nb(ct.shp), style="B")
# rownames(adjCT) <- colnames(adjCT) <- as.character(ct.shp$ZCTA5CE10)
# adjCT <- adjCT[order(rownames(adjCT)), order(rownames(adjCT))]
# save(adjCT, file="/Users/aritrahalder/Desktop/new_sp_tw_1/adjCT.RData")
load("adjCT.RData")

idCT <- match(ct_zip$zip,rownames(adjCT)) 
ct_zip <- ct_zip[!is.na(idCT),]
degCT <- diag(rowSums(adjCT)); rownames(degCT) <- colnames(degCT) <- rownames(adjCT)

# Simulation parameters
N <- c(1e4,2e4,3e4,5e4)
psiL=c(7,11,24,100)
psiU=c(11,24,100,200)

# smooth
splits <- seq(min(ct_zip$latitude), max(ct_zip$latitude), length.out = 10)
splits_val <- seq(-3, 3, length.out = 9)
speff <- sapply(ct_zip$latitude, function(x){
  if(x<=splits[2]) return(-3)
  for(i in 2:8) if(x<=splits[(i+1)] & x>splits[i]) return(splits_val[i])
  if(x>=splits[9]) return(3)
})

col_vec <- RColorBrewer::brewer.pal(9,"GnBu")
col <- sapply(speff, function(x){
  for(i in 1:9) if(x == splits_val[i]) return(col_vec[i])
})
speff <- speff-mean(speff)

for(i in 1:length(N)){
  for(j in 1:length(psiL)){
    set.seed (1234)
    mean_cov  <- data.frame(cbind(rbinom(N[i],1,0.5),
                                  rbinom(N[i],4,0.5),
                                  rnorm(N[i],0,1),
                                  rnorm(N[i],0,1)))
    disp_cov <- data.frame(cbind(rbinom(N[i],1,0.5),
                                 rbinom(N[i],4,0.5),
                                 rnorm(N[i],0,0.1),
                                 rnorm(N[i],0,0.1)))
    colnames(mean_cov)  <- paste("X",1:4,sep=""); colnames(disp_cov) <-  paste("Z",1:4,sep="")
    mean_cov[,1:2] <- apply(mean_cov[,1:2],2,function(x) as.factor(x))
    disp_cov[,1:2] <- apply(disp_cov[,1:2],2,function(x) as.factor(x))
    
    x_mat <- cbind(as.numeric(mean_cov[,1]=="1"),
                   as.numeric(mean_cov[,2]=="1"),
                   as.numeric(mean_cov[,2]=="2"),
                   as.numeric(mean_cov[,2]=="3"),
                   as.numeric(mean_cov[,2]=="4"),
                   mean_cov[,3],
                   mean_cov[,4])
    x_mat <- cbind(1,x_mat)
    z_mat <- cbind(as.numeric(mean_cov[,1]=="1"),
                   as.numeric(mean_cov[,2]=="1"),
                   as.numeric(mean_cov[,2]=="2"),
                   as.numeric(mean_cov[,2]=="3"),
                   as.numeric(mean_cov[,2]=="4"),
                   mean_cov[,3],
                   mean_cov[,4])
    z_mat <- cbind(1,z_mat)
    set.seed(NULL)
    
    # Spatial Pattern:: Blockwise
    # N: 10000, 30000, 50000
    # prop.zero: 0.15, 0.30, 0.60, 0.80
    
    
    
    nsamp <- rep(floor(N[i]/nrow(ct_zip)),nrow(ct_zip))
    nsamp[sample(1:nrow(ct_zip),size = (N[i]-nrow(ct_zip)*floor(N[i]/nrow(ct_zip))))] <- nsamp[sample(1:nrow(ct_zip),size = (N[i]-nrow(ct_zip)*floor(N[i]/nrow(ct_zip))))]+1
    
    thetaHat <- rnorm(N[i],-0.16,0.02)
    theta <- thetaHat*exp(-rep(speff,nsamp)/2)
    muHat <- 4/thetaHat^2
    mu_sim <- 4/theta^2
    phi_sim <- runif(N[i],psiL[j],psiU[j])
    gammaTrue <- crossprod(solve(crossprod(z_mat,z_mat)),crossprod(z_mat,log(phi_sim)))
    betaTrue <- crossprod(solve(crossprod(x_mat,x_mat)),crossprod(x_mat,log(muHat)))
    mean_cov$zip <- rep(rownames(adjCT),nsamp)
    mean_cov$speff <- rep(speff,nsamp)
    
    # simulated response
    y <- rtweedie(N[i],xi = 1.5,mu = mu_sim,phi = phi_sim)
    
    # final data
    sim_data <- data.frame(longitude = rep(ct_zip$longitude,nsamp),
                           latitude = rep(ct_zip$latitude,nsamp),
                           mean_cov,
                           mu_sim,
                           phi_sim,
                           y)
    print(round(length(y[y==0])/N[i],2))
    # boxplot(y~round(sim_data$speff,2), xlab="Spatial Effect")
    obj <- list(data=sim_data,
                betaT=betaTrue,
                gammaT=gammaTrue,
                alphaT=speff,
                X=x_mat,
                Z=z_mat)
    save(obj, file=paste("/Users/aritrahalder/Desktop/paper1_code/sim_data/smooth_",i,j,".RData",sep=""))
  }
}
