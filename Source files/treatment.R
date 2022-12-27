### Prepare Datasets for prediction ###
# Inputs:
# - y = response variable
# - X = regressors matrix
# - col.m: column index of metric variables, e.g. col.m=c(1:5)
# - col.n: column index of non-metric variables, e.g. col.n=c(6:10)
# Output:
# - Databases with treated non-metric variables following four methods (NMC, MNMC,
# One-hot Encoding and original data set with no treatment).
# Scripts to be used:
# - "normal_mean_coding.R": normal mean coding method for non-metric variables treatment
source("data_processing/normal_mean_coding.R")
# - "mnmc.R": Multivariate normal mean coding method for non-metric variables treatment
require("combinat")
require("nleqslv")
require("calibrate")
source("data_processing/mnmc.R")
source("data_processing/auxiliaryfunMNMC.R") # auxiliary functions of "mnmc.R"
# - "onehot.R": One hot Encoding method for non-metric variables treatment
source("data_processing/onehot.R")
source("data_processing/accesories.R")
source("data_processing/predictmnmc.R")
library("plyr")

## FUNCTIONS:

preprocessing <- function(Xtrain,Xtest, col.m, col.n){
  X.m<-Xtrain[, col.m] # metric regressors
  X.n<-Xtrain[, col.n] # non-metric regressors
  
  X.m_test<-Xtest[, col.m] # metric regressors
  X.n_test<-Xtest[, col.n] 
  
  n<-dim(Xtrain)[1]; k<-dim(Xtrain)[2] # data dimensions
  
  #para volver a acomodar
  npm=c(col.n,col.m)
  names(npm)<-seq(npm)
  npm=sort(npm)
  new_order=as.numeric(names(npm))

  ### Non-metric treatment methodologies
  # 1. (univariate) Normal Mean Coding
  # 2. Multavariate Normal Mean Coding
  # 3. One-hot encoding

  # 1. (univariate) Normal Mean Coding
    # Treatment on train set
    temp<-normal.mean.coding(X=X.n)
    Xnmc_train_aux <-cbind(temp, X.m)
    Xnmc_train <-Xnmc_train_aux[,new_order]
    # Treatment on test set
    Xnmc_test_nm <- cbind()
    for (col in 1:length(col.n)){
      values_ordinal <- c(sort(unique(X.n[,col])))
      values_nmc <- c(sort(unique(temp[,col])))
      mapped_values <- mapvalues(X.n_test[,col], from = values_ordinal, to = values_nmc)
      Xnmc_test_nm <- cbind(Xnmc_test_nm,mapped_values)
    }
    colnames(Xnmc_test_nm) <- NULL
    Xnmc_test <- cbind(Xnmc_test_nm,X.m_test)
    NMC <- list(Xnmc_train,Xnmc_test)
    
  
  # 2. Multivariate Normal Mean Coding
    # Treatment on train set
    temp <- X.n+1
    Delta0=NULL
    temp<-mnmc(X=temp, W=X.m, Delta0 = Delta0, ConvCriteria=2,tol=.01)
    Xmnmc_train_aux <-temp[[1]] # treated train set
    p = length(col.n) # number of non-metric variables
    Xmnmc_train=Xmnmc_train_aux[,new_order]
    Delta_aux = temp[[2]]
    Delta=Delta_aux[new_order,new_order]
    theta = temp[[3]]
    
    
    # Treatment on test set
    Vtest = cbind(X.n_test+1,X.m_test) 
    Xmnmc_test_aux=predictmnmc(Vtest,p,Delta,theta)
    Xmnmc_test=Xmnmc_test_aux[,new_order]
    MNMC = list(Xmnmc_train,Xmnmc_test)
  
    
    
    
  
  # 3. One-hot encoding. OJO estos no los devulve en el orden original
    # temporal rbind() of train and test to do the onehot encoding
    r_train = dim(X.n)[1] # rows in train
    r_test = dim(X.n_test)[1] # rows in test
    temp = rbind(X.n,X.n_test)
    colnames(temp)<-1:dim(temp)[2]
    temp = one_hot_encoding(as.data.frame(temp), columns = names(as.data.frame(temp)))
    temp_train = temp[(1:r_train),]
    temp_test = temp[((r_train+1):dim(temp)[1]),]
    Xohe_train<-cbind(temp_train,X.m)
    Xohe_test <- cbind(temp_test,X.m_test)
    OHE = list(Xohe_train,Xohe_test)
    
  # 4. No treatment
    Xnt_train <- Xtrain[,new_order]
    Xnt_test <- Xtest[,new_order]
    NT <- list(Xnt_train,Xnt_test)
  
  
  # Return all preprocessed variables
  return(list(NMC=NMC,MNMC=MNMC,OneHot=OHE,NoTreat=NT))
}

preprocessing_binary <- function(Xtrain,Xtest, col.m, col.n, col.b){
  X.m<-Xtrain[, col.m] # metric regressors
  X.n<-Xtrain[, col.n] # non-metric regressors
  X.b <- Xtrain[, col.b] # binary regressors
  X.m_test<-Xtest[, col.m] # metric regressors
  X.n_test<-Xtest[, col.n] # non-metric regressors
  X.b_test <- Xtest[,col.b] # binary regressors
  
  n<-dim(Xtrain)[1]; k<-dim(Xtrain)[2] # data dimensions
  
  #para volver a acomodar
  npm=c(col.n,col.m)
  names(npm)<-seq(npm)
  npm=sort(npm)
  new_order=as.numeric(names(npm))
  
  ### Non-metric treatment methodologies
  # 1. (univariate) Normal Mean Coding
  # 2. Multavariate Normal Mean Coding
  # 3. One-hot encoding
  
  # 1. (univariate) Normal Mean Coding
  # Treatment on train set
  temp<-normal.mean.coding(X=X.n)
  Xnmc_train_aux <-cbind(temp, X.m)
  Xnmc_train <-Xnmc_train_aux[,new_order]
  Xnmc_train <- cbind(Xnmc_train,X.b)
  # Treatment on test set
  Xnmc_test_nm <- cbind()
  for (col in 1:length(col.n)){
    values_ordinal <- c(sort(unique(X.n[,col])))
    values_nmc <- c(sort(unique(temp[,col])))
    mapped_values <- mapvalues(X.n_test[,col], from = values_ordinal, to = values_nmc)
    Xnmc_test_nm <- cbind(Xnmc_test_nm,mapped_values)
  }
  colnames(Xnmc_test_nm) <- NULL
  Xnmc_test <- cbind(Xnmc_test_nm,X.m_test,X.b_test)
  NMC <- list(Xnmc_train,Xnmc_test)
  
  
  # 2. Multivariate Normal Mean Coding
  # Treatment on train set
  temp <- X.n+1
  Delta0=NULL
  temp<-mnmc(X=temp, W=X.m, Delta0 = Delta0, ConvCriteria=2,tol=.01)
  Xmnmc_train_aux <-temp[[1]] # treated train set
  p = length(col.n) # number of non-metric variables
  Xmnmc_train=Xmnmc_train_aux[,new_order]
  Xmnmc_train = cbind(Xmnmc_train,X.b)
  Delta_aux = temp[[2]]
  Delta=Delta_aux[new_order,new_order]
  theta = temp[[3]]
  
  
  # Treatment on test set
  Vtest = cbind(X.n_test+1,X.m_test) 
  Xmnmc_test_aux=predictmnmc(Vtest,p,Delta,theta)
  Xmnmc_test=Xmnmc_test_aux[,new_order]
  Xmnmc_test=cbind(Xmnmc_test,X.b_test)
  MNMC = list(Xmnmc_train,Xmnmc_test)
  
  # 3. One-hot encoding. OJO estos no los devulve en el orden original
  # temporal rbind() of train and test to do the onehot encoding
  r_train = dim(X.n)[1] # rows in train
  r_test = dim(X.n_test)[1] # rows in test
  temp = rbind(X.n,X.n_test)
  colnames(temp)<-1:dim(temp)[2]
  temp = one_hot_encoding(as.data.frame(temp), columns = names(as.data.frame(temp)))
  temp_train = temp[(1:r_train),]
  temp_test = temp[((r_train+1):dim(temp)[1]),]
  Xohe_train<-cbind(temp_train,X.m,X.b)
  Xohe_test <- cbind(temp_test,X.m_test,X.b_test)
  OHE = list(Xohe_train,Xohe_test)
  
  # 4. No treatment
  Xnt_train <- Xtrain[,new_order]
  Xnt_train <- cbind(Xnt_train,X.b)
  Xnt_test <- Xtest[,new_order]
  Xnt_test <- cbind(Xnt_test,X.b_test)
  NT <- list(Xnt_train,Xnt_test)
  
  
  # Return all preprocessed variables
  return(list(NMC=NMC,MNMC=MNMC,OneHot=OHE,NoTreat=NT))
}

preprocessing_binary <- function(Xtrain,Xtest, col.m, col.n, col.b){
  X.m<-Xtrain[, col.m] # metric regressors
  X.n<-Xtrain[, col.n] # non-metric regressors
  X.b <- Xtrain[, col.b] # binary regressors
  X.m_test<-Xtest[, col.m] # metric regressors
  X.n_test<-Xtest[, col.n] # non-metric regressors
  X.b_test <- Xtest[,col.b] # binary regressors
  
  n<-dim(Xtrain)[1]; k<-dim(Xtrain)[2] # data dimensions
  
  #para volver a acomodar
  npm=c(col.n,col.m)
  names(npm)<-seq(npm)
  npm=sort(npm)
  new_order=as.numeric(names(npm))
  
  ### Non-metric treatment methodologies
  # 1. (univariate) Normal Mean Coding
  # 2. Multavariate Normal Mean Coding
  # 3. One-hot encoding
  
  # 1. (univariate) Normal Mean Coding
  # Treatment on train set
  temp<-normal.mean.coding(X=X.n)
  Xnmc_train_aux <-cbind(temp, X.m)
  Xnmc_train <-Xnmc_train_aux[,new_order]
  Xnmc_train <- cbind(Xnmc_train,X.b)
  # Treatment on test set
  Xnmc_test_nm <- cbind()
  for (col in 1:length(col.n)){
    values_ordinal <- c(sort(unique(X.n[,col])))
    values_nmc <- c(sort(unique(temp[,col])))
    mapped_values <- mapvalues(X.n_test[,col], from = values_ordinal, to = values_nmc)
    Xnmc_test_nm <- cbind(Xnmc_test_nm,mapped_values)
  }
  colnames(Xnmc_test_nm) <- NULL
  Xnmc_test <- cbind(Xnmc_test_nm,X.m_test,X.b_test)
  NMC <- list(Xnmc_train,Xnmc_test)
  
  
  # 2. Multivariate Normal Mean Coding
  # Treatment on train set
  temp <- X.n+1
  Delta0=NULL
  temp<-mnmc(X=temp, W=X.m, Delta0 = Delta0, ConvCriteria=2,tol=.01)
  Xmnmc_train_aux <-temp[[1]] # treated train set
  p = length(col.n) # number of non-metric variables
  Xmnmc_train=Xmnmc_train_aux[,new_order]
  Xmnmc_train = cbind(Xmnmc_train,X.b)
  Delta_aux = temp[[2]]
  Delta=Delta_aux[new_order,new_order]
  theta = temp[[3]]
  
  
  # Treatment on test set
  Vtest = cbind(X.n_test+1,X.m_test) 
  Xmnmc_test_aux=predictmnmc(Vtest,p,Delta,theta)
  Xmnmc_test=Xmnmc_test_aux[,new_order]
  Xmnmc_test=cbind(Xmnmc_test,X.b_test)
  MNMC = list(Xmnmc_train,Xmnmc_test)
  
  # 3. One-hot encoding. OJO estos no los devulve en el orden original
  # temporal rbind() of train and test to do the onehot encoding
  r_train = dim(X.n)[1] # rows in train
  r_test = dim(X.n_test)[1] # rows in test
  temp = rbind(X.n,X.n_test)
  colnames(temp)<-1:dim(temp)[2]
  temp = one_hot_encoding(as.data.frame(temp), columns = names(as.data.frame(temp)))
  temp_train = temp[(1:r_train),]
  temp_test = temp[((r_train+1):dim(temp)[1]),]
  Xohe_train<-cbind(temp_train,X.m,X.b)
  Xohe_test <- cbind(temp_test,X.m_test,X.b_test)
  OHE = list(Xohe_train,Xohe_test)
  
  # 4. No treatment
  Xnt_train <- Xtrain[,new_order]
  Xnt_train <- cbind(Xnt_train,X.b)
  Xnt_test <- Xtest[,new_order]
  Xnt_test <- cbind(Xnt_test,X.b_test)
  NT <- list(Xnt_train,Xnt_test)
  
  
  # Return all preprocessed variables
  return(list(NMC=NMC,MNMC=MNMC,OneHot=OHE,NoTreat=NT))
}



test_treatment <- function(treatment,Xtrain,Xtest, binary_treatment){
  if(binary_treatment == TRUE){
    X.m<-Xtrain[, col.m] # metric regressors
    X.n<-Xtrain[, col.n] # non-metric regressors
  
    X.m_test<-Xtest[, col.m] # metric regressors
    X.n_test<-Xtest[, col.n] 
  
    n<-dim(Xtrain)[1]; k<-dim(Xtrain)[2] # data dimensions
  
    #para volver a acomodar
    npm=c(col.n,col.m)
    names(npm)<-seq(npm)
    npm=sort(npm)
    new_order=as.numeric(names(npm))
  
    ### Non-metric treatment methodologies
    # 1. (univariate) Normal Mean Coding
    # 2. Multavariate Normal Mean Coding
    # 3. One-hot encoding
    if(treatment == "NM"){
      # 1. (univariate) Normal Mean Coding
      # Treatment on train set
      temp<-normal.mean.coding(X=X.n)
      Xnmc_train_aux <-cbind(temp, X.m)
      Xnmc_train <-Xnmc_train_aux[,new_order]
      # Treatment on test set
      Xnmc_test_nm <- cbind()
      for (col in 1:length(col.n)){
        values_ordinal <- c(sort(unique(X.n[,col])))
        values_nmc <- c(sort(unique(temp[,col])))
        mapped_values <- mapvalues(X.n_test[,col], from = values_ordinal, to = values_nmc)
        Xnmc_test_nm <- cbind(Xnmc_test_nm,mapped_values)
    }
      colnames(Xnmc_test_nm) <- NULL
      Xnmc_test <- cbind(Xnmc_test_nm,X.m_test)
      results <- list(Xnmc_train,Xnmc_test)
  }
    if(treatment == "MNMC"){
      # 2. Multivariate Normal Mean Coding
      # Treatment on train set
      temp <- X.n+1
      Delta0=NULL
      temp<-mnmc(X=temp, W=X.m, Delta0 = Delta0, ConvCriteria=2,tol=.01)
      Xmnmc_train_aux <-temp[[1]] # treated train set
      p = length(col.n) # number of non-metric variables
      Xmnmc_train=Xmnmc_train_aux[,new_order]
      Delta_aux = temp[[2]]
      Delta=Delta_aux[new_order,new_order]
      theta = temp[[3]]
    
    
      # Treatment on test set
      Vtest = cbind(X.n_test+1,X.m_test) 
      Xmnmc_test_aux=predictmnmc(Vtest,p,Delta,theta)
      Xmnmc_test=Xmnmc_test_aux[,new_order]
      results = list(Xmnmc_train,Xmnmc_test)
  }
  
  
  
    if(treatment == "OHE"){
      # 3. One-hot encoding. OJO estos no los devulve en el orden original
      # temporal rbind() of train and test to do the onehot encoding
      r_train = dim(X.n)[1] # rows in train
      r_test = dim(X.n_test)[1] # rows in test
      temp = rbind(X.n,X.n_test)
      colnames(temp)<-1:dim(temp)[2]
      temp = one_hot_encoding(as.data.frame(temp), columns = names(as.data.frame(temp)))
      temp_train = temp[(1:r_train),]
      temp_test = temp[((r_train+1):dim(temp)[1]),]
      Xohe_train<-cbind(temp_train,X.m)
      Xohe_test <- cbind(temp_test,X.m_test)
      results = list(Xohe_train,Xohe_test)
  }
  
  
    if (treatment == "NT cont" | treatment == "NT mixed"){
      # 4. No treatment
      Xnt_train <- Xtrain[,new_order]
      Xnt_test <- Xtest[,new_order]
      results <- list(Xnt_train,Xnt_test)
  }
  
  # Return all preprocessed variables
  return(results)
  }else{
    X.m<-Xtrain[, col.m] # metric regressors
    X.n<-Xtrain[, col.n] # non-metric regressors
    X.b <- Xtrain[, col.b] # binary regressors
    X.m_test<-Xtest[, col.m] # metric regressors
    X.n_test<-Xtest[, col.n] # non-metric regressors
    X.b_test <- Xtest[,col.b] # binary regressors
    
    n<-dim(Xtrain)[1]; k<-dim(Xtrain)[2] # data dimensions
    
    #para volver a acomodar
    npm=c(col.n,col.m)
    names(npm)<-seq(npm)
    npm=sort(npm)
    new_order=as.numeric(names(npm))
    
    ### Non-metric treatment methodologies
    # 1. (univariate) Normal Mean Coding
    # 2. Multavariate Normal Mean Coding
    # 3. One-hot encoding
    if(treatment == "NM"){
      # 1. (univariate) Normal Mean Coding
      # Treatment on train set
      temp<-normal.mean.coding(X=X.n)
      Xnmc_train_aux <-cbind(temp, X.m)
      Xnmc_train <-Xnmc_train_aux[,new_order]
      Xnmc_train <- cbind(Xnmc_train,X.b)
      # Treatment on test set
      Xnmc_test_nm <- cbind()
      for (col in 1:length(col.n)){
        values_ordinal <- c(sort(unique(X.n[,col])))
        values_nmc <- c(sort(unique(temp[,col])))
        mapped_values <- mapvalues(X.n_test[,col], from = values_ordinal, to = values_nmc)
        Xnmc_test_nm <- cbind(Xnmc_test_nm,mapped_values)
      }
      colnames(Xnmc_test_nm) <- NULL
      Xnmc_test <- cbind(Xnmc_test_nm,X.m_test,X.b_test)
      results <- list(Xnmc_train,Xnmc_test)
    }
    
    if(treatment == "MNMC"){
      # 2. Multivariate Normal Mean Coding
      # Treatment on train set
      temp <- X.n+1
      Delta0=NULL
      temp<-mnmc(X=temp, W=X.m, Delta0 = Delta0, ConvCriteria=2,tol=.01)
      Xmnmc_train_aux <-temp[[1]] # treated train set
      p = length(col.n) # number of non-metric variables
      Xmnmc_train=Xmnmc_train_aux[,new_order]
      Xmnmc_train = cbind(Xmnmc_train,X.b)
      Delta_aux = temp[[2]]
      Delta=Delta_aux[new_order,new_order]
      theta = temp[[3]]
    
    
      # Treatment on test set
      Vtest = cbind(X.n_test+1,X.m_test) 
      Xmnmc_test_aux=predictmnmc(Vtest,p,Delta,theta)
      Xmnmc_test=Xmnmc_test_aux[,new_order]
      Xmnmc_test=cbind(Xmnmc_test,X.b_test)
      results = list(Xmnmc_train,Xmnmc_test)
    }
    
    if(treatment == "OHE"){
      # 3. One-hot encoding. OJO estos no los devulve en el orden original
      # temporal rbind() of train and test to do the onehot encoding
      r_train = dim(X.n)[1] # rows in train
      r_test = dim(X.n_test)[1] # rows in test
      temp = rbind(X.n,X.n_test)
      colnames(temp)<-1:dim(temp)[2]
      temp = one_hot_encoding(as.data.frame(temp), columns = names(as.data.frame(temp)))
      temp_train = temp[(1:r_train),]
      temp_test = temp[((r_train+1):dim(temp)[1]),]
      Xohe_train<-cbind(temp_train,X.m,X.b)
      Xohe_test <- cbind(temp_test,X.m_test,X.b_test)
      results = list(Xohe_train,Xohe_test)
    }
    if (treatment == "NT cont" | treatment == "NT mixed"){
      # 4. No treatment
      Xnt_train <- Xtrain[,new_order]
      Xnt_train <- cbind(Xnt_train,X.b)
      Xnt_test <- Xtest[,new_order]
      Xnt_test <- cbind(Xnt_test,X.b_test)
      results <- list(Xnt_train,Xnt_test)
    }
    
    # Return all preprocessed variables
    return(results)
  }
}

