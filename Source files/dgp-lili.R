# simulation ###############################################################

# programs for discretization #################################################

### discretize###########################
# This program discretizes a vector to a discrete variables 
# with certain number of unique values.

# inputs 
# x.star: metric vector
# m.j: number of unique values

# outputs
# x: discretized variable

discretize<-function(x.star, m.j){
  # define objects
  n<-length(x.star)
  
  # By very low chance the thresholds may overlap or be at the lowest or highest values. 
  # It leads to the number of unique values lower than desired and hence I avoid it via loop.
  repeat{
    # thresholds
    tau<-sort(runif(n=(m.j-1), min = 0, max = 1))
    
    # the values of x corresponding to the thresholds in the quantile of the empirical CDF
    tau.star<-sort(x.star)[ceiling(tau*n)]
    
    # discretize x.star  
    x<-matrix(unlist(lapply(as.list(tau.star), function(a) x.star>a)), ncol = m.j-1, byrow = FALSE)
    x<-apply(x, 1, sum)
    
    # break upon the desired number of unique categories
    if(length(unique(x))==m.j){break}
  }
  
  # report
  return(x)
}

### discretize.matrix ##################
# This program discretizes a matrix to a discrete variables 
# with certain number of unique values.

# inputs 
# X.star: metric matrix
# m: number of unique values

# outputs
# X: discretized variable

discretize.matrix<-function(X.star, m){
  # define objects
  n<-dim(X.star)[1]; k<-dim(X.star)[2]
  m<-as.list(m)
  X.star<-split(X.star, rep(1:ncol(X.star), each = nrow(X.star)))
  
  # discretize
  wl<-as.list(1:k)
  X<-matrix(unlist(lapply(wl, function(i) discretize(x.star=X.star[[i]], m.j=m[[i]]))), ncol = k, byrow = FALSE)
  
  # report
  return(X)
}

# end: programs for discretization


# programs for data generating processes (DGP)########################################

### gm1 ###############################################
# This program generates metric variables following the first and third DGP.

# inputs
# n=100: number of observations
# k=10: number of variables
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.delta=1: sigma^2 parameter for the distribution of delta

# outputs
# X.star: metric variables
# xi.1: latent variable

gm1<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.delta=1){
  # define objects
  lambda<-as.matrix(rep(sqrt(1/k), k))
  
  # generate xi.1
  if(dist.xi==1){
    xi.1<-as.matrix(rnorm(n=n, mean=par.xi.1[1], sd=sqrt(par.xi.1[2])))
  } else {
    xi.1<-as.matrix(rlnorm(n=n, meanlog=par.xi.1[1], sdlog=sqrt(par.xi.1[2])))
  }
  
  # generate Delta
  delta<-matrix(rnorm(n=n*k, mean=0, sd=sqrt(par.delta)), nrow=n, ncol=k)
  
  # generate X.star
  X.star<-xi.1%*%t(lambda)+delta
  
  # report
  return(list(X.star=X.star, xi.1=xi.1))
}

### gm2 ###################################
# This program generates metric variables following the second and fourth DGP.

# inputs
# n=100: number of observations
# k=10: number of variables (It has to be a multiple of 5, e.g. 5, 10, 15,...)
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.xi.2=c(0,5): mu and sigma^2 parameter for the distribution of xi.2
# par.delta=1: sigma^2 parameter for the distribution of delta

# outputs
# X.star: metric variables
# xi.1: latent variable

gm2<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1){
  # define objects
  lambda.1<-as.matrix(rep(sqrt(1/k), k)) 
  lambda.2<-as.matrix(rep(c(1/(2*sqrt(k)),1/(2*sqrt(k)),1/(2*sqrt(k)),1/(2*sqrt(k)),-2/(sqrt(k))), k/5))
  
  # generate xi.1 and xi.2
  if(dist.xi==1){
    xi.1<-as.matrix(rnorm(n=n, mean=par.xi.1[1], sd=sqrt(par.xi.1[2])))
    xi.2<-as.matrix(rnorm(n=n, mean=par.xi.2[1], sd=sqrt(par.xi.2[2])))
  } else {
    xi.1<-as.matrix(rlnorm(n=n, meanlog=par.xi.1[1], sdlog=sqrt(par.xi.1[2])))
    xi.2<-as.matrix(rlnorm(n=n, meanlog=par.xi.2[1], sdlog=sqrt(par.xi.2[2])))
  }
  
  # generate Delta
  delta<-matrix(rnorm(n=n*k, mean=0, sd=sqrt(par.delta)), nrow=n, ncol=k)
  
  # generate X.star
  
  X.star<-3.5*xi.1%*%t(lambda.1)+0.3*xi.2%*%t(lambda.2)+delta
  
  # report
  return(list(X.star=X.star, xi.1=xi.1, xi.2=xi.2))
}

# con los scores no ortogonales
gm2no<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1){
  # define objects1
  lambda.1<-as.matrix(rep(sqrt(1/k), k)) 
  lambda.2<-as.matrix(c(rep(1/(sqrt(k/2)), floor(k/2)), rep(8, k-floor(k/2)))*4) 

  # generate xi.1 and xi.2
  if(dist.xi==1){
    xi.1<-as.matrix(rnorm(n=n, mean=par.xi.1[1], sd=sqrt(par.xi.1[2])))
    xi.2<-as.matrix(rnorm(n=n, mean=par.xi.2[1], sd=sqrt(par.xi.2[2])))
  } else {
    xi.1<-as.matrix(rlnorm(n=n, meanlog=par.xi.1[1], sdlog=sqrt(par.xi.1[2])))
    xi.2<-as.matrix(rlnorm(n=n, meanlog=par.xi.2[1], sdlog=sqrt(par.xi.2[2])))
  }
  
  # generate Delta
  delta<-matrix(rnorm(n=n*k, mean=0, sd=sqrt(par.delta)), nrow=n, ncol=k)
  
  # generate X.star
  X.star<-3.5*xi.1%*%t(lambda.1)+0.3*xi.2%*%t(lambda.2)+delta
  
  # report
  return(list(X.star=X.star, xi.1=xi.1, xi.2=xi.2))
}

### dgp1: 1 factor latente relacionado con y y relacionado con x ###############
# This program generates data following the first DGP.

# inputs
# n=100: number of observations
# k=10: number of variables
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.delta=1: sigma^2 parameter for the distribution of delta
# par.e=1: sigma^2 parameter for the distribution of E

# outputs
# X.star: metric variables
# y: outcome variable
# Xb: true fit

dgp1<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.delta=1, par.e=0.01){
  # generate latent variables and metric variables
  yeon<-gm1(n=n, k=k, dist.xi=dist.xi, par.xi.1=par.xi.1, par.delta=par.delta/(9*k))

  # latent variable
  xi.1<-yeon$xi.1
  # xi.1<-unit.variance(yeon$xi.1) # debug
  
  # metric variables
  X.star<-yeon$X.star
  
  # betas
  betas<-as.matrix(1)
  
  # error term
  e<-rnorm(n=n, mean=0, sd=sqrt(par.e))
  
  # generate outcome variables
  Xb<-xi.1%*%betas
  #Xb<-unit.variance(Xb) # debug
  y<-Xb+e
  
  # report
  result<-list(y=y, X.star=X.star, Xb=Xb)  
  return(result)
}

### dgp2: 1 factor latente relacionado con y con x, y otro factor relacionado con x ortogonal con el primero ######
# This program generates data following the second DGP.

# inputs
# n=100: number of observations
# k=10: number of variables (It has to be a multiple of 5, e.g. 5, 10, 15,...)
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.xi.2=c(0,5): mu and sigma^2 parameter for the distribution of xi.2
# par.delta=1: sigma^2 parameter for the distribution of delta
# par.e=1: sigma^2 parameter for the distribution of E

# outputs
# X.star: metric variables
# y: outcome variable
# Xb: true fit

dgp2<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1, par.e=0.01){
  # generate latent variables and metric variables
  yeon<-gm2(n=n, k=k, dist.xi=dist.xi, par.xi.1=par.xi.1, par.xi.2=par.xi.2, par.delta=par.delta/(9*k))
  
  # latent variable
  xi.1<-yeon$xi.1
  # xi.1<-unit.variance(yeon$xi.1) # debug
  
  # metric variables  
  X.star<-yeon$X.star
  
  # betas
  betas<-as.matrix(1)
  
  # error term
  e<-rnorm(n=n, mean=0, sd=sqrt(par.e))
  
  # generate outcome variables
  Xb<-xi.1%*%betas
 # Xb<-unit.variance(Xb) # debug
  y<-Xb+e
  
  # report
  result<-list(y=y, X.star=X.star, Xb=Xb)  
  return(result)
}

### dgp3: 2 factores latentes, ambos relacionados con x y con y , y son ortogonales ######
# This program generates data following the third DGP.

# inputs
# n=100: number of observations
# k=10: number of variables
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.delta=1: sigma^2 parameter for the distribution of delta
# par.e=1: sigma^2 parameter for the distribution of E

# outputs
# X.star: metric variables
# y: outcome variable
# Xb: true fit

dgp3<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1, par.e=0.01){
  # generate latent variables and metric variables
  yeon<-gm2(n=n, k=k, dist.xi=dist.xi, par.xi.1=par.xi.1, par.xi.2=par.xi.2, par.delta=par.delta/(9*k))
  
  # latent variable
  xi.1<-yeon$xi.1
  xi.2<-yeon$xi.2

  # xi.1<-unit.variance(yeon$xi.1) # debug
  
  # metric variables  
  X.star<-yeon$X.star
  
  # betas
 
  
  # error term
  e<-rnorm(n=n, mean=0, sd=sqrt(par.e))
  betas=as.matrix(c(1,1))
  # generate outcome variables
  Xb<-cbind(xi.1,xi.2)%*%betas
  
  # Xb<-unit.variance(Xb) # debug
  y<-Xb+e
  
  # report
  result<-list(y=y, X.star=X.star, Xb=Xb)  
  return(result)
}

### dgp4: 1 factor latente relacionado con x e y, y un factor relacionado sólo con x, no ortogonal con el primero ###
# This program generates data following the fourth DGP.

# inputs
# n=100: number of observations
# k=10: number of variables (It has to be a multiple of 5, e.g. 5, 10, 15,...)
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.xi.2=c(0,5): mu and sigma^2 parameter for the distribution of xi.2
# par.delta=1: sigma^2 parameter for the distribution of delta
# par.e=1: sigma^2 parameter for the distribution of E

# outputs
# X.star: metric variables
# y: outcome variable
# Xb: true fit


#hacemos que y dependa de par.xi2 y par.xi1
dgp4<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1, par.e=0.01){
  # generate latent variables and metric variables
  yeon<-gm2no(n=n, k=k, dist.xi=dist.xi, par.xi.1=par.xi.1, par.xi.2=par.xi.2, par.delta=par.delta/(9*k))
  
  # latent variable
  xi.1<-yeon$xi.1
  # xi.1<-unit.variance(yeon$xi.1) # debug
  
  # metric variables  
  X.star<-yeon$X.star
  
  # betas
  betas<-as.matrix(1)
  
  # error term
  e<-rnorm(n=n, mean=0, sd=sqrt(par.e))
  
  # generate outcome variables
  Xb<-xi.1%*%betas
  # Xb<-unit.variance(Xb) # debug
  y<-Xb+e
  
  # report
  result<-list(y=y, X.star=X.star, Xb=Xb)  
  return(result)
}

### dgp5: x e y relacionadas con 2 factores latentes no ortogonales entre sí ########
# This program generates data following the third DGP.

# inputs
# n=100: number of observations
# k=10: number of variables
# dist.xi=1: If dist.xi=1, xi is normal distributed. 
# Otherwise, xi is log normal distributed.
# par.xi.1=c(0,1): mu and sigma^2 parameter for the distribution of xi.1
# par.delta=1: sigma^2 parameter for the distribution of delta
# par.e=1: sigma^2 parameter for the distribution of E

# outputs
# X.star: metric variables
# y: outcome variable
# Xb: true fit

dgp5<-function(n=100, k=10, dist.xi=1, par.xi.1=c(0,1), par.xi.2=c(0,1), par.delta=1, par.e=0.01){
  # generate latent variables and metric variables
  yeon<-gm2no(n=n, k=k, dist.xi=dist.xi, par.xi.1=par.xi.1, par.xi.2=par.xi.2, par.delta=par.delta/(9*k))
  
  # latent variable
  xi.1<-yeon$xi.1
  xi.2<-yeon$xi.2
  # xi.1<-unit.variance(yeon$xi.1) # debug
  
  # metric variables  
  X.star<-yeon$X.star
  
  # betas
  
  
  # error term
  e<-rnorm(n=n, mean=0, sd=sqrt(par.e))
  betas=as.matrix(rep(1,2))
  # generate outcome variables
  Xb<-cbind(xi.1,xi.2)%*%betas
  #Xb<-X.star%*%betas
  # Xb<-unit.variance(Xb) # debug
  y<-Xb+e
  
  # report
  result<-list(y=y, X.star=X.star, Xb=Xb)  
  return(result)
}




