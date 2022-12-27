#This file contains auxialliary functions used in the main function mnmc 
#(multivariate normal mean coding)


#### 1. findThresholds ####

library(pracma) #I use this function to uses same functions as Matlab code.
#May be we can delete it.
library(rootSolve) #Check if it is required


findThresholds<-function(X,DeltaX){
  # This function estimate the thresholds of p ordinal variables (X)
  #assuming an underlying normal distribution for p latent variables (Z)
  # X are p ordinal response variables (n \times p)
  # DeltaX is the covariance matrix of Z|X (p \times p)
  
  X=t(t(X)) #To read X as matrix
  n=dim(X)[1] #Number of observations
  p=dim(X)[2] #Number of ordinal variables
  unos=rep(1,n) #vector of ones used in objective function
  
  Theta=matrix(list(), p, 1) #structure to store the thresholds
  for (j in 1:p){ 
    medias=rep(0,p) 
    if (p==1){dj<-DeltaX
    Xj<-X} #Case of only one ordinal variable
    else
    {dj<-DeltaX[j,j] 
    Xj<-X[,j]}
    Kj<-max(Xj)
    theta_j = matrix(0,nrow=1,ncol=Kj-1)
    
    for (k in 1:(Kj-1)){
      njk<-length(which(Xj<=k)) #count observations with category values 
      #less or equal to the k-th category
      
      objfun <- function(x) {
        temp <- (njk - sum(pnorm((x*unos)/dj))); 
        return(temp);
      } #Objective function
      
      x0=c(min(medias)-5*dj,max(medias)+5*dj); #range to search roots
      x0<-na.omit(x0)
      if (length(x0) != 2) {x0=c(-5,5)}
      
     
temp=uniroot(objfun,x0,extendInt="yes", maxiter=500) #we use uniroot function to
# find the thresholds (roots) of objfun
      
      theta_j[k]=temp$root
      
    }
    Theta[[j]]=theta_j
  }
  return(Theta);
}


##### 2. updateEzANDEzz #####

updateEzANDEzz<-function(Ez,Ezz,Theta,X,Delta,idx){
  i = idx[1]; j= idx[2];
  X=t(t(X));
  Xij=X[i,j];
  thetaj = Theta[[j]];
  kj = length(thetaj)+1;
  if (is.null(dim(Delta))){sigmatilde_ij = sqrt(Delta)
  } else { 
    Sigma_jmj = Delta[j,]; Sigma_jmj=Sigma_jmj[-j];
  Sigma_mjmj = Delta; Sigma_mjmj=Sigma_mjmj[-j,-j]; 
  Sigma_jmj=as.matrix(Sigma_jmj)
  sigmatilde_ij = sqrt(Delta[j,j] - t(Sigma_jmj)%*%solve(Sigma_mjmj)%*%(Sigma_jmj));
  }
  if (is.null(dim(Delta))){Emutilde_ij=0 
  } else {
  Ez_imj = Ez[i,]; Ez_imj = Ez_imj[-j];
  Ez_imj=as.matrix(Ez_imj)
  Emutilde_ij = t(Sigma_jmj)%*%solve(Sigma_mjmj)%*%(Ez_imj);
  Ezz_imjmj = cbind(Ezz[i,,1:size(Ezz,3)]); Ezz_imjmj=Ezz_imjmj[-j,-j];
  }
  
  if (Xij==1){
    deltatilde_ij = (thetaj[Xij] - Emutilde_ij)/sigmatilde_ij;
    coefic = -1*dnorm(deltatilde_ij)/pnorm(deltatilde_ij);
    if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
    #Updating E(Zij|all)
    Eznew_ij = Emutilde_ij + coefic * sigmatilde_ij;
    
    #Updating E(Zij*Zij | all)     
    coef2 = (0 - deltatilde_ij*dnorm(deltatilde_ij))/(pnorm(deltatilde_ij) - 0);
    if (is.infinite(coef2)||is.nan(coef2)){coef2=0}
    
    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 +  2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;
    
  } else if (Xij==kj){
    deltatilde_ijm = (thetaj[Xij-1] - Emutilde_ij)/sigmatilde_ij;
    coefic = dnorm(deltatilde_ijm)/(1-pnorm(deltatilde_ijm));
    if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
    Eznew_ij = Emutilde_ij + coefic * sigmatilde_ij;
    
    
    coef2 = (deltatilde_ijm*dnorm(deltatilde_ijm) - 0)/(1 - pnorm(deltatilde_ijm));
    if (is.infinite(coef2)||is.nan(coef2)){coef2=0}
    
    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 +  2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;
    
  } else { deltatilde_ijm = (thetaj[Xij-1] - Emutilde_ij)/sigmatilde_ij;
  deltatilde_ij = (thetaj[Xij] - Emutilde_ij)/sigmatilde_ij;
  coefic = (dnorm(deltatilde_ijm) - dnorm(deltatilde_ij)) / (pnorm(deltatilde_ij) - pnorm(deltatilde_ijm));
  if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
  #Update E(Zij|all)
  Eznew_ij = Emutilde_ij + coefic*sigmatilde_ij;
  
  #Update E(Zij*Zij | all)         
  coef2 = (deltatilde_ijm*dnorm(deltatilde_ijm) - deltatilde_ij*dnorm(deltatilde_ij))/(pnorm(deltatilde_ij) - pnorm(deltatilde_ijm));
  if (is.infinite(coef2)||is.nan(coef2)){coef2=0}
  
  Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 + 2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;
  
  }
  Output=list(Eznew_ij,Ezznew_ijj);
  return(Output);
}


##### 3. computeEzANDEzz #####

computeEzANDEzz<-function(V,p,Delta,Theta){
  V=t(t(V))
  t=size(V,2)
  n=size(V,1)
  X = V[,1:p];
  if (t==p){
    Ez = zeros(n,t);
    Ezz = array(0,c(n,t,t)); 
  } else {
    X = V[,1:p];
    W = V[,(p+1):t]
    Ez = zeros(n,t); 
    Ez[,(p+1):t] = W;
    #Ezz = array(0,c(n,p,p));
   Ezz = array(0,c(n,t,t)); 
    for (j in 1:n) {
      Ezz[j,(p+1):t,(p+1):t]=cov(W)}
    
 }
  
  SS=zeros(p,p)
 # SS=zeros(p,p)
  Eznew = t(t(Ez));
  Ezznew =(Ezz);
  
  # Initialization
  StopNotMet1 = 1;
  StopNotMet2 = 1;
  
  history = 1;
  Sold = eye(p);
  iter = 0;
  
  while(StopNotMet1==1 && StopNotMet2==1 && iter<100){
    iter = iter + 1;
    for (i in 1:n){
      for (j in 1:p){
        Auxx = updateEzANDEzz(Ez,Ezz,Theta,V,Delta,c(i,j));
        a1=as.numeric(Auxx[1])
        a2=as.numeric(Auxx[2])
        Eznew[i,j]= a1
        Ezznew[i,j,j]=a2
      }
      for (j in 1:p){
        for (jprima in 1:p){
          if (j!=jprima){
            Ezznew[i,j,jprima] = Eznew[i,j]*Eznew[i,jprima]}}
      }
    }
    # Update M and S
    
    for (j in 1:p){
      for (jprima in 1:p){
        Ezznew[,j,jprima]=as.vector(Ezznew[,j,jprima])
        SS[j,jprima] = sum(Ezznew[,j,jprima])}}
    
    M = Eznew;
    
    Ezz = Ezznew;
    Ez = Eznew;
    

    outConv1 = checkConvergenceS(SS,Sold,tol=NULL);
    StopNotMet1=outConv1
    outConv2 = checkConvergenceEz(Eznew,Ez,tol=NULL);
    StopNotMet2=outConv2
    #history=outConv[[2]]
    Sold = SS;
  }  
  resultado=list(M,SS)
  return(resultado)
}

#####Functions used for convergence####

##### 4. loglik #####

loglik<-function(DeltaOld,Delta,n){
  # This function compute log-likelihood when updating Delta
  invDelta = inv(Delta)
 valfun= -(1/2)*n* log(det(Delta))-n/2*sum(diag((invDelta%*%DeltaOld)))
 return(valfun)
}

##### 5. checkConvergence #####

checkConvergence<- function(DeltaOld, Delta,n, history,tol){
if (is.null(tol)){tol=1e-10} 
StopNotMet = 1;
funval = loglik(DeltaOld, Delta,n);
if(is.na(funval)){StopNotMet=0} else if ((funval < history)||( abs(history-funval)/abs(funval) < tol)){
StopNotMet = 0};
history = funval;
out<-list(StopNotMet,funval)
return(out)
}


##### 6. checkConvergenceS #####

checkConvergenceS<-function(S,Sold,tol){
  if (is.null(tol)){tol=1e-5} 
  StopNotMet = 1;
  if(abs(norm(S-Sold, type='F')/norm(Sold, type='F'))<tol){StopNotMet=0} else{StopNotMet=1}
  return(StopNotMet)
}


##### 7. checkConvergenceEz #####

checkConvergenceEz<-function(Eznew,Ez,tol){
  if (is.null(tol)){tol=1e-5} 
  StopNotMet = 1;
  if(abs(norm(Eznew-Ez, type='F')/norm(Ez, type='F'))<tol){StopNotMet=0} else{StopNotMet=1}
  return(StopNotMet)
}

