###############################################################
# normal mean coding
###############################################################

# Normal means: Kolenikov and Angeles [2009] discuss a rescaling of ordinal variable. See p137

# normal.mean.coding.vector
# input
# y: ordinal variable

# output
# result: normal mean coded variable
# quantification

normal.mean.coding.vector<-function(y){
  #y<-round(y, 5)
  
  # alpha_hat: It returns the threshold of an ordinal variable x
  # input
  # x: ordinal variable
  
  # output
  # result: threshold
  
  alpha_hat<-function(x){
    j<-unique(x)
    j<-j[order(j)]
    N<-length(x)
    
    result<-c(0)
    for(a in 1:length(j)){
      n<-length(x[x==j[a]|x<j[a]])
      result<-append(result, qnorm((-0.5+n)/N, mean = 0, sd = 1, log = FALSE))
    }
    result<-result[-1]
    
    return(result)
  }
  
  alpha_hat_y<-alpha_hat(y)
  wpnm<-dnorm(append(-Inf, alpha_hat_y[-length(alpha_hat_y)]))-dnorm(alpha_hat_y)
  wpnm<-cbind(unique(y)[order(unique(y))], wpnm)
  
  result<-c(0) 
  for(i in 1:length(y)){ 
    result<-append(result, vlookup(y[i], wpnm))
  }
  result<-result[-1]
  
  return(list(result=result, quantification=wpnm))
}

# normal.mean.coding
# input
# X: ordinal variable data

# output
# X: normal mean coded ordinal variable data

normal.mean.coding<-function(X){
  X<-as.matrix(X)
  
  p<-dim(X)[2]
  for(j in 1:p){
    X[,j]<-normal.mean.coding.vector(X[,j])$result
  }
  return(X)
}

# normal.mean.coding.prediction
# input
# X_train: ordinal variable data from training set
# X_test: ordinal variable data from testing set

# output
# X_train: quantified ordinal variable data from training set
# X_test: quantified ordinal variable data from testing set

normal.mean.coding.prediction<-function(X_train, X_test){
  X_train<-as.matrix(X_train); X_test<-as.matrix(X_test)
  
  X_test<-round(X_test, 5)
  p_train<-dim(X_train)[2]
  
  for(j in 1:p_train){
    temp<-normal.mean.coding.vector(X_train[,j])
    X_train[,j]<-temp$result
    
    wpnm<-temp$quantification
    result<-c(0) 
    for(i in 1:length(X_test[,j])){ 
      result<-append(result, vlookup(X_test[i,j], wpnm))
    }
    result<-result[-1]
    
    X_test[,j]<-result
  }
  return(list(X_train=X_train, X_test=X_test))
}
