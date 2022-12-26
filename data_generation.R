### Data Generation ###
# Inputs:
# - n = number of observations
# - k = number of variables
# - pnm = proportion of non-metric variables


# Scripts to be used:
#   Data generation
# - "dgp-lili.R": generates data
source("data_generation/dgp-lili.R")
#   Data treatment
# - "tm.R": returns the number of categories of the non-metric variables
source("data_generation/tm.R")
library(foreign)


## FUNCTIONS:

call_dgp <- function(n,k,pnm = 0.5, dgp){
  nnmv <- round(k*pnm, 0) # number of nonmetric variables
  col.m <- 1:(k-nnmv) # metric column index
  col.n <- (k-nnmv+1):k # non-metric column index
  enc = 7
  m<-tm(enc=enc, nnmv=nnmv) # number of unique categories (referenced outside)
  if(dgp==1){
    temp <- dgp1(n=n, k=k) # default values for other arguments in "dgp-lili.R"
  }else if(dgp==2){
    temp <- dgp2(n=n, k=k)
  }else if(dgp==3){
    temp <- dgp3(n=n, k=k)
  }else if(dgp==4){
    temp <- dgp4(n=n, k=k)
  }else if(dgp==5){
    temp <- dgp5(n=n, k=k)
  }else{temp <- dgp6(n=n, k=k)}
  # Results of dgp
  X <- temp$X.star # Latent continuos
  X[,col.n] <- discretize.matrix(X.star=temp$X.star[,col.n], m=m) # discretize non-metric variables
  y <- temp$y # regressand
  Xb <- temp$Xb # true fitted values
  # unit variance scaling for metric variables  
  inv.sigma <- sqrt(1/apply(X[,col.m], 2, var))
  diag.inv.sigma <- matrix(0, ncol=length(inv.sigma), nrow=length(inv.sigma))
  diag(diag.inv.sigma) <- inv.sigma
  X[,col.m]<- X[,col.m] %*% diag.inv.sigma
  
  # Output
  return(list(y=y, X=X, col.m=col.m, col.n=col.n))
}


call_dgp_binary <- function(n,k,pnm = 0.5, pbin=0.3, dgp){
  nnmv <- round(k*pnm, 0) # number of nonmetric variables
  nbin <- round(k*pbin, 0)
  col.m <- 1:(k-nnmv-nbin) # metric column index
  col.n <- (k-nnmv-nbin+1):(k-nbin) # non-metric column index
  col.b <- (k-nbin+1):k # binary column index
  enc = 7
  m<-tm(enc=enc, nnmv=nnmv) # number of unique categories (referenced outside)
  mbin<-tm_binary(enc=2, nbin=nbin) 
  if(dgp==1){
    temp <- dgp1(n=n, k=k) # default values for other arguments in "dgp-lili.R"
  }else if(dgp==2){
    temp <- dgp2(n=n, k=k)
  }else if(dgp==3){
    temp <- dgp3(n=n, k=k)
  }else if(dgp==4){
    temp <- dgp4(n=n, k=k)
  }else if(dgp==5){
    temp <- dgp5(n=n, k=k)
  }else{temp <- dgp6(n=n, k=k)}
  # Results of dgp
  X <- temp$X.star # Latent continuos
  X[,col.n] <- discretize.matrix(X.star=temp$X.star[,col.n], m=m) # discretize non-metric variables
  X[,col.b] <- discretize.matrix(X.star=temp$X.star[,col.b], m=mbin) # discretize non-metric variables
  y <- temp$y # regressand
  Xb <- temp$Xb # true fitted values
  # unit variance scaling for metric variables  
  inv.sigma <- sqrt(1/apply(X[,col.m], 2, var))
  diag.inv.sigma <- matrix(0, ncol=length(inv.sigma), nrow=length(inv.sigma))
  diag(diag.inv.sigma) <- inv.sigma
  X[,col.m]<- X[,col.m] %*% diag.inv.sigma
  
  # Output
  return(list(y=y, X=X, col.m=col.m, col.n=col.n, col.b=col.b))
}

