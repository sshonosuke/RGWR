###----------------------------------------------------###
###          Robust GWR with gamma-divergence          ###
###----------------------------------------------------###


### Function for fitting robust GWR with fixed tuning parameters  
## Input
# Y: vector of response variables
# X: (n,p)-matrix of covariates (should not include an intercept)
# Sp: (n,2)-matrix of location information (longitude and latitude)
# gam: robustness parameter (default: gam=0.3)
# b: spatial bandwidth (default: b=0.1)
# st.error: If 'T', standard errors are calculated.
# maxit: maximum number of iterations of MM-algorithm (default: maxit=100)
## Output (list object)
# beta: (n,p+1)-matrix of estimated location-wise regression coefficients 
# sig: vector of estimated location-wise error variances
# SE: (n,p+1)-matrix of standard errors of coefficients 

RGWR <- function(Y, X, Sp, gam=0.3, b=0.1, st.error=T, maxit=100){
  n <- length(Y)
  p <- dim(X)[2]
  XX <- cbind(1, X)
  dd <- as.matrix( dist(Sp) )
  W <- exp(-0.5*dd^2/b^2)
  hBeta <- matrix(NA, n, p+1)
  hSig <- c()
  for(i in 1:n){
    est <- MM(y=Y, x=X, ww=W[i,], gam=gam, maxit=maxit)
    hBeta[i,] <- est$beta
    hSig[i] <- est$sig
  }
  SE <- matrix(NA, n, p+1)
  if(st.error){
    for(i in 1:n){
      ww <- W[i,]
      beta <- hBeta[i,]
      mu <- as.vector(XX%*%beta)
      resid <- Y - mu
      sig <- hSig[i]
      uu <- dnorm(Y, mu, sig)^(gam)
      EE <- ww*uu*resid*XX
      I <- t(XX)%*%(ww^2*uu^2*resid^2*XX)
      resid2 <- gam*resid^2/sig^2 - 1
      IJ <- solve( t(XX)%*%(ww*uu*resid2*XX) )
      SE[i,] <- sqrt( diag(IJ%*%I%*%IJ) )
    }
  }
  return( list(beta=hBeta, sig=hSig, SE=SE) )
}






### Function for selecting two tuning parameters in RGWR  
## Input
# Y: vector of response variables
# X: (n,p)-matrix of covariates (should not include an intercept)
# Sp: (n,2)-matrix of location information (longitude and latitude)
# b.max: maximum value of spatial bandwidth 
# L: the number of candidate values for bandwidth 
# gam.set: candidate set for gamma (robustness parameter)
## Output (list object)
# gam: selected robustness parameter
# b: selected bandwidth
# HS: candidate values of gamma and their H-scores 
# CV: candidate values of bandwidth and their CV-values

RGWR.sel <- function(Y, X, Sp, b.max=NULL, L=20, gam.set=NULL){
  ## candidate sets
  if(is.null(b.max)){ 
    b.max <- 0.5*median(dist(Sp)) 
  }
  b.set <- seq(0, b.max, length=L+1)[-1]
  if(is.null(gam.set)){ 
    gam.set <- c(0, 0.01, 0.03, seq(0.05, 0.5, by=0.05))
  }
  J <- length(gam.set)
  
  ## selection of gamma
  Score <- c()
  for(j in 1:J){
    gam <- gam.set[j]
    Fit <- RGWR(Y, X, Sp, gam=gam, b=b.max)
    hBeta <- Fit$beta
    hMu <- apply(cbind(1,X)*hBeta, 1, sum)
    hSig <- Fit$sig
    ww <- dnorm(Y, hMu, hSig)^gam
    C <- ( sqrt(1+gam)*(2*pi*hSig^2)^(gam/2) )^(gam/(1+gam))
    H <- 2*C*(gam*(Y-hMu)^2-hSig^2)*ww/hSig^4 + C^2*ww^2*(Y-hMu)^2/hSig^4 
    Score[j] <- sum(H)
  }
  hgam <- gam.set[which.min(Score)]
  
  ## selection of bandwidth
  val <- c()
  for(l in 1:L){
    val[l] <- RCV(Y, X, Sp, b=b.set[l], gam=hgam)
  }
  b.opt <- b.set[which.max(val)]
  
  ## Summary 
  Result <- list(gam=hgam, b=b.opt, HS=cbind(gam.set, Score), CV=cbind(b.set, val)) 
  return( Result )
}







###    Function for MM-algorithm     ###
MM <- function(y, x, ww, gam=0.3, maxit=100){
  ## initial value 
  nn <- length(y)
  p <- dim(x)[2]
  if(sum(ww)>p+1){
    init <- lm(y~x, weights=ww)
    hBeta <- coef(init)
    hSig <- summary(init)$sigma
  }else{
    init <- lm(y~x)
    hBeta <- coef(init)
    hSig <- summary(init)$sigma
  }
  ## iteration
  XX <- cbind(1, x)
  for(k in 1:maxit){
    # weight
    Mu <- as.vector(XX%*%hBeta)
    dens <- dnorm(y, Mu, hSig)^(gam)
    uu <- ww*dens
    U <- uu/sum(uu)   # normalized weight 
    # Beta
    if( sum(U>1/nn)>(p+2) ){
      new.Beta <- as.vector( solve(t(XX)%*%(U*XX))%*%t(XX)%*%(U*y) )
    }else{
      new.Beta <- hBeta
    }
    new.Mu <- as.vector(XX%*%new.Beta)
    # sigma
    hSig <- sqrt( (1+gam)*sum(U*(y-new.Mu)^2) )
    # convergence check
    dd <- 100*sum(abs(new.Beta-hBeta))/sum(abs(hBeta)+0.0001)
    hBeta <- new.Beta
    if(dd<0.01){ break() }
  }
  ## Summary
  return( list(beta=hBeta, sig=hSig) )
}



###    Function for robust cross-validation    ###
RCV <- function(Y, X, Sp, gam=0.3, b=0.1, maxit=10){
  ## preparation
  n <- length(Y)
  p <- dim(X)[2]
  XX <- cbind(1, X)
  dd <- as.matrix( dist(Sp) )
  W <- exp(-0.5*dd^2/b^2)
  ## iteration
  hMu <- c()
  hSig <- c()
  for(i in 1:n){
    est <- MM(y=Y[-i], x=X[-i,], ww=W[i,-i], gam=gam, maxit=maxit)
    hMu[i] <- sum(XX[i,]*est$beta)
    hSig[i] <- est$sig
  }
  ## computation of CV
  if(gam>0){
    wdens <- dnorm(Y, hMu, hSig)^(gam)
    val <- log(mean(wdens))/gam + 0.5*gam/(1+gam)*log(mean(hSig^2))
    val <- n*exp(gam*val)/gam
  }
  if(gam==0){
    val <- sum(dnorm(Y, hMu, hSig, log=T))
  }
  return(val)
}


