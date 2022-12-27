###############################################################
# Accessories
# centralize: a program to centralize data
# unit.length: rescale columns as unit length
# unit.variance: rescale columns as unit variance
# vlookup: a program analogous to vlookup in Microsoft Excel
# cv.msep: it calculates 10 fold cross validation mean squared error.
# svd2: it calculates the first singular vectors and value with 
# computationally expensive but stable way.
# interaction: It builds interaction terms between clustering variables
# and regressors.
###############################################################
###############################################################
library(foreign)
library(mgcv)
#library(quantreg)
library(nnet)
library(combinat)
library(xtable)
library(calibrate)
library(vegan)
library(pracma)
library(mvtnorm)
library(nleqslv)
# library(MatrixModels): This package is not easy to be installed on the server and I do not remember the use.
# Delete it later if there is no problem.
library(polycor) # a package to calculate polychoric correlation

###############################################################
# centralize
###############################################################

# a program to centralize data

# input
# target: data

# output
# target: centralized data

centralize<-function(target){
  target<-as.matrix(target)
  means<-colMeans(target)
  iota<-rep(1, dim(target)[1])
  target<-target-iota%*%t(means)
  return(target)
}

###############################################################
# end: centralize
###############################################################

###############################################################
# unit.length
###############################################################

# a program to make the columns of a matrix of unit length 

# input
# target: data

# output
# target: data with each column recaled as unit length

ss<-function(x){sum(x^2, na.rm=TRUE)}# ss returns sum of squares of a vector

unit.length<-function(target){
  target<-as.matrix(target)
  inv.sigma<-sqrt(1/apply(target, 2, ss))
  diag.inv.sigma<-matrix(0, ncol=length(inv.sigma), nrow=length(inv.sigma))
  diag(diag.inv.sigma)<-inv.sigma
  target<-target%*%diag.inv.sigma
  return(target)
}

###############################################################
# end: unit.length
###############################################################

###############################################################
# unit.variance
###############################################################

# a program to make the columns of a matrix of unit variance

# input
# target: data

# output
# target: data with each column recaled as unit length

unit.variance<-function(target){
  target<-as.matrix(target)
  inv.sigma<-sqrt(1/apply(target, 2, var))
  diag.inv.sigma<-matrix(0, ncol=length(inv.sigma), nrow=length(inv.sigma))
  diag(diag.inv.sigma)<-inv.sigma
  target<-target%*%diag.inv.sigma
  return(target)
}

###############################################################
# end: unit.variance
###############################################################

###############################################################
# vlookup
###############################################################

# a program analogous to vlookup in Microsoft Excel

# input
# val: the value to be matched
# mat: the vlookup matrix

# output
# result: the matched value

vlookup<-function(val, mat){
  if(sum(mat[,1] == val)!=0){
    return(mat[mat[,1] == val, 2])
  } else {
    return(NA)
  }
}

###############################################################
# end: vlookup
###############################################################

###############################################################
# cv.msep
###############################################################

# cv.msep calculates 10 fold cross validation mean squared error.
# inputs
# y: regressand
# X: regressors
# logit_link=FALSE : logit_link=TRUE triggers logit model instead of linear model
# id=NA : One can specify segments in cross-validataion. This is useful when
# one wants to create bootstrap confidence interval. id should be a vector
# of size being equal to the size of y and include interger values from 1 to 10.

# outputs
# msep: MSEP

cv.msep<-function(y, X, logit_link=FALSE, id=NA){
  # all inputs as matrix
  y<-as.matrix(y)
  X<-as.matrix(X)
  
  # define objects
  N<-length(y)
  if(sum(!is.na(id))!=0){}else{id <- sample(1:10, dim(X)[1], replace=TRUE)} 
  
  # calculate msep
  msep<-0
  for(k in 1:10){
    # Define inputs
    wy=y[id!=k]
    wX=X[id!=k,]
    wy0=y[id==k]
    wX0=X[id==k,]
    
    # calculate fitted values
    if(logit_link==FALSE){
      mod<-lm(wy~wX)
      iota0<-rep(1, dim(as.matrix(wX0))[1])
      b<-mod$coefficients
      if(sum(is.na(b))>0){print(paste("Warning: Regressors are not of full rank at the ", k, "th segment",". Unidentified coefficients are replaced by zero.", sep=""))}
      b[is.na(b)]<-0 # Substitute NA coefficients due to multicollinearity to zero.
      fk<-cbind(iota0, wX0)%*%b
    } else {
      workingdata<-data.frame(wy, wX)
      workingdata0<-data.frame(wy0, wX0)
      colnames(workingdata0)<-colnames(workingdata)
      mod<-glm(formula = wy ~ wX, family = binomial(logit), data = workingdata)
      fk<-predict(mod, newdata = workingdata0,  type = "response")
    }
    
    # msep for kth segment
    Nk<-length(fk)    
    msep<-msep+(sum((fk-y[id==k])^2, na.rm=T)/Nk)/10
  }
  
  # report result
  return(msep)
}

###############################################################
# end: cv.msep
###############################################################

###############################################################
# f.parameters.xi
###############################################################

# This program calculates the parameters of log normal distributed xi
# as a function of variance and skewness.

# inputs
# v=1: variance
# s=13: skewness

# outputs
# x: x=c(log mean parameter, log variance parameter)

f.parameters.xi<-function(v=1, s=13){
  # skekur
  # This function returns c(variance-v, skewness-s) of a log normal distribution
  # as a function of x1 and x2.
  
  # input
  # x: a vector containing c(x1, x2)
  
  # output 
  # y: a vector containing c(variance-v, skewness-s)
  
  skekur <- function(x) {
    y <- numeric(2)
    y[1]<-(exp(x[2])-1)*exp(2*x[1]+x[2])-v
    y[2]<-(exp(x[2])+2)*sqrt(exp(x[2])-1)-s
    y
  }
  
  xstart <- c(0, 1) # arbitrary starting points
  
  x<-nleqslv(xstart, skekur, control=list(btol=.01))[[1]]
  
  # report
  return(x)
}

###############################################################
# end: f.parameters.xi
###############################################################

###############################################################
# f.parameters.delta
###############################################################

# This program calculates the variance of delta as a function of signal noise ratio

# inputs
# alpha=3: signal noise ratio
# k: number of variables

# outputs
# result: variance of delta.i

f.parameters.delta<-function(alpha=3, k) 1/((alpha^2)*k)

###############################################################
# end: f.parameters.delta
###############################################################

###############################################################
# svd2
###############################################################

# svd sometimes fails to deliver results. I guess that they use an algorithm
# with cheaper computational costs and less stability.
# I write down an stable algorithm, but maybe with more computational costs.

# inputs
# X: data matrix

# outputs
# d: a vector containing the first singular values of x
# u: a matrix whose columns contain the first left singular vectors of x
# v: a matrix whose columns contain the first right singular vectors of x

svd2<-function(X){
  temp<-eigen(t(X)%*%X)
  d<-sqrt(temp$values[1])
  v<-as.matrix(temp$vectors[,1])
  u<-as.matrix(unit.length(X%*%v))
  
  return(list(d=d, u=u, v=v))
}

###############################################################
# end: svd2
###############################################################

###############################################################
# interaction
###############################################################

# This program builds interaction terms between clustering variables
# and regressors.

# input
# X: regressors
# G: clustering variables (in dummy coding)

# output
# result: a matrix containing regressors, intercepts corresponding to the
# clustering variables and interactions between the regressors and the regressors.

interaction<-function(X, G){
  # transform objects as matrices
  X<-as.matrix(X)
  G<-as.matrix(G)
  
  # add column names in case there is none present
  if(is.null(colnames(X))){colnames(X)<-paste("X", 1:dim(X)[2], sep="")}
  if(is.null(colnames(G))){colnames(G)<-paste("G", 1:dim(G)[2], sep="")}
  
  # an object to be reported
  result<-cbind(X, G)
  
  # make interaction terms and add to the result
  for(j in 1:dim(G)[2]){
    temp<-X*G[,j]
    colnames(temp)<-paste(colnames(X), ".", colnames(G)[j])
    result<-cbind(result, temp)
  }
  
  # report
  return(result)
}

###############################################################
# end: interaction
###############################################################
