#This function estimate the conditional expectation of ordinal variables X,
#i.e. M=E(Z|X). Here we follow the first step algorithm of Forzani et.al
#(2018).
#INPUTS:
# X: ordinal predictor matrix, each row is an observation and each column a predictor variable.
# W: continuous predictor matrix, each row is an observation and each column a predictor variable.
#ConvCriteria =1 if we take distance


# OUTPUT:
# M: first moment estimates for the latent variables, given the observed Parameter Estimates: Delta, Theta; data.

source("data_processing/auxiliaryfunMNMC.R")

mnmc<-function(X,W=NULL,Delta0=NULL, ConvCriteria =2, tol=0.01){
    
 
   if (is.null(X)){M=W;  Delta=cov(W); Theta=NULL} else{
  
    X=t(t(X)) ; # To ensure that X is taken as a matrix
    p = dim(X)[2]; # Number of ordinal variables
    n = dim(X)[1]; # Number of Observations
    
    
    if(is.null(W)){dimen=p} else{dimen=p+dim(W)[2]} #Total number of variables
    if(is.null(Delta0)){Delta0 =diag(rep(1,dimen))} # Initial Delta (Delta0)
    #Delta0 = Identity as Default

  if(is.null(W)){V=X} else { Wo = scale(W,center=TRUE,scale=FALSE);
V=cbind(X,Wo)} #Centering continuous variables

  #Initialization
  Delta=Delta0;
  history = -1e6; #history in likelihood values
  StopNotMet = 1; # Zero means aright to the stop criteria, and one otherwise
  iter = 0; # For iterations

  while (StopNotMet==1 && iter<20){
  iter = iter+1;
  print(iter)
  
  # Update thresholds
  Theta = findThresholds(X,Delta);
  DeltaOld=Delta;

  # Ez and Ezz computation
  MSOut=computeEzANDEzz(V,p,Delta,Theta)

    Maux=MSOut[[1]] #Estimated M=[E(Z|X) Wo]
  SSaux= MSOut[[2]] #Estimated S=E(Z%*%t(Z)|X)
  Maux=t(t(Maux)) #To read M as matrix
  SSaux=t(t(SSaux)) # To read S as matrix
  MauxO = Maux[,1:p]; #Here we take only the estimated latent normal var E(Z|X)
  
  # We construct the covariance matrix 
  if (is.null(W)){
    Delta1 = (SSaux)*(1/n);
    Delta = .5*(Delta1+t(Delta1)); # To ensure symmetry of Delta
    
    VV=MauxO
    f=diag(c((diag(Delta[1:p,1:p]))^(-1/2))) 
    Delta = f%*%Delta%*%f #Normalizing the Covariance matrix to have correlation 
    #for ordinal block (needed for identification)  
    M=Maux
    } else {
  Delta1 = (SSaux)*(1/n);
  cross=1/n*t(MauxO)%*%Wo;
  first=cbind(Delta1,cross)
  second= cbind(t(cross),cov(W))
  Delta= rbind(first,second) 
  Delta = .5*(Delta+t(Delta)); # To ensure symmetry of Delta
  
  VV=cbind(MauxO,W)
  f=diag(c((diag(Delta[1:p,1:p]))^(-1/2),rep(1, dim(W)[2]))) 
  Delta = f%*%Delta%*%f #Normalizing the Covariance matrix to have correlation 
  #for ordinal block (needed for identification)
  }  
        
  
    
  
  #Checking convergence with log-likelihood
  
  if (ConvCriteria ==1){
  checkOut = checkConvergence(DeltaOld, Delta,n,history,tol)} else {
  checkOut=checkConvergenceS(DeltaOld,Delta,tol)}

  
  
 StopNotMet = unlist(checkOut[1])
 history = unlist(checkOut[2])
 
 
 M = cbind(Maux[,1:p],W);   
 
  }
 
  }
  result=list(M,Delta, Theta)
  
  return(result)
    
}

