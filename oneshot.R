rm(list=ls())

##  simulation  (one-shot)    ##
set.seed(2)

library(spgwr)   
library(GWmodel)
library(MASS)
source("RGWR-function.R")


##  settings
n <- 500    # number of samples
om <- 0.03   # contamination ratio 


# generation of sampling locations
Sp <- matrix(NA, n, 2)
for(i in 1:n){
  dd <- 0
  while(dd<(0.5^2)){
    rn <- c(runif(1,-1,1),runif(1,0,2))
    dd <- sum(rn[1]^2+0.5*rn[2]^2)
  }
  Sp[i,] <- rn
}
plot(Sp)


#  covariates
phi <- 0.4   # range parameter for covariates
dd <- as.matrix(dist(Sp))
mat <- exp(-dd/phi)

z1 <- mvrnorm(1, rep(0, n), mat)
z2 <- mvrnorm(1, rep(0, n), mat)
x1 <- z1
rr <- 0.75
x2 <- rr*z1 + sqrt(1-rr^2)*z2
X <- cbind(x1, x2)


# Spatially varying parameters 
V1 <- exp(-dd/1)
V2 <- exp(-dd/2)
V3 <- exp(-dd/3)
Beta0 <- mvrnorm(1, rep(0,n), 2*V1)
Beta1 <- mvrnorm(1, rep(0,n), 2*V2)
Beta2 <- mvrnorm(1, rep(0,n), 2*V3)
Sig.true <- 1
Beta.true <- cbind(Beta0, Beta1, Beta2)


## Data
Mu <- apply(cbind(1,X)*Beta.true, 1, sum)
Y <- rnorm(n, Mu, Sig.true)
ch <- rbinom(n, 1, om)
Y[ch==1] <- Y[ch==1] + 5
plot(Y)



## standard GWR
b.opt <- gwr.sel(Y~X, coords=Sp, verbose=F)
b.opt
fit <- gwr(Y~X, coords=Sp, bandwidth=b.opt, se.fit=T, hatmatrix=T)
est.GWR <- cbind(fit$SDF$`(Intercept)`, fit$SDF$`Xx1`, fit$SDF$`Xx2`)


## proposed robust GWR
opt <- RGWR.sel(Y, X, Sp)
opt$gam
opt$b
rfit <- RGWR(Y, X, Sp, gam=opt$gam, b=opt$b)
est.RGWR <- rfit$beta


# MSE
apply((est.GWR - Beta.true)^2, 2, mean) 
apply((est.RGWR - Beta.true)^2, 2, mean) 



# plot
par(mfcol=c(1, 2))
for(k in 2:3){
  Est <- cbind(Beta.true[,k], est.GWR[,k], est.RGWR[,k])
  ran <- range(Est)
  plot(Est[,1:2], ylim=ran, xlim=ran, ylab="estimate", 
       xlab="true", main=paste0("x", k-1), pch=16)
  points(Est[,c(1,3)], col=2, pch=16)
  abline(0, 1)
  legend("topleft", legend=c("GWR", "RGWR"), col=1:2, pch=16)
}
par(mfcol=c(1, 1))



