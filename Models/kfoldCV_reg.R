### k-fold CV for regression ###
# Inputs:
# - dim_reduction (integer): number of directions PCA and PLS
# - number_of_folds (integer): number of K-folds
# - binary_treatment (logical): if binary variables are included or not in non-metric variables treatment
kfoldCV_reg <- function(dim_reduction=1,number_of_folds=5,binary_treatment=FALSE,ids){
  K <- number_of_folds
  a <- dim_reduction
  ids <- id
  # Create empty array to be filled by results:
  # 1st dimension: prediction methods (3 - linear model, non-parametric, inverse regression)
  # 2nd dimension: treatment and reduction methods (10)
  # 3rd dimension: number of folds
  Kresults <- array(NA, dim=c(K,10,3))
  for(i in 1:K){
    flag <- paste0("Running the ",i,"-th fold",sep="")
    print(flag)
    ytrain <- y[id!=i]
    ytest <- y[id==i]
    Xtrain <- X[id!=i,]
    Xtest <- X[id==i,]
    if(binary_treatment){
      #Projections with treatment of binary variables
      projections <- pca.pls(ytrain,ytest,Xtrain,Xtest,col.m,col.n,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }else{#Projections with NO treatment of binary variables
      projections <- pca.pls_binary(ytrain,ytest,Xtrain,Xtest,col.m,col.n,col.b,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }
    for (pr in 1:length(proj_test)) {
      # Linear Regression
      fit <- linear_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MSE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 1] <- MSE_test
      # Non-Parametric Regression
      fit <- np_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 2] <- MSE_test
      # Inverse Regression
      fit <- inverse_reg(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 3] <- MSE_test
    }
  }
  results = colMeans(Kresults, dim=1)
  colnames(results) <- c("Linear Reg","Non-Parametric Reg", "Inverse Reg")
  rownames(results) <- c("PCA NM","PCA MNMC", "PCA OHE", "PCA NT cont", "PCA NT mixed",
                         "PLS NM","PLS MNMC", "PLS OHE", "PLS NT cont", "PLS NT mixed")
  return(results)
}



kfoldCV_reg2 <- function(dim_reduction=1,number_of_folds=5,binary_treatment=FALSE,ids){
  K <- number_of_folds
  a <- dim_reduction
  ids <- id
  # Create empty array to be filled by results:
  # 1st dimension: prediction methods (3 - linear model, non-parametric, inverse regression)
  # 2nd dimension: treatment and reduction methods (10)
  # 3rd dimension: number of folds
  Kresults <- array(NA, dim=c(K,10,1))
  for(i in 1:K){
    flag <- paste0("Running the ",i,"-th fold",sep="")
    print(flag)
    ytrain <- y[id!=i]
    ytest <- y[id==i]
    Xtrain <- X[id!=i,]
    Xtest <- X[id==i,]
    if(binary_treatment){
      #Projections with treatment of binary variables
      projections <- pca.pls(ytrain,ytest,Xtrain,Xtest,col.m,col.n,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }else{#Projections with NO treatment of binary variables
      projections <- pca.pls_binary(ytrain,ytest,Xtrain,Xtest,col.m,col.n,col.b,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }
    for (pr in 1:length(proj_test)) {
      # Linear Regression
      fit <- linear_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MSE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 1] <- MSE_test
      # Non-Parametric Regression
      #fit <- np_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      #MSE_test <-  fit$MSE_TEST
      #Kresults[i , pr , 2] <- MSE_test
      # Inverse Regression
      #fit <- inverse_reg(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      #MSE_test <-  fit$MSE_TEST
      #Kresults[i , pr , 3] <- MSE_test
    }
  }
  results = colMeans(Kresults, dim=1)
  colnames(results) <- c("Linear Reg")#,"Non-Parametric Reg")#, "Inverse Reg")
  rownames(results) <- c("PCA NM","PCA MNMC", "PCA OHE", "PCA NT cont", "PCA NT mixed",
                         "PLS NM","PLS MNMC", "PLS OHE", "PLS NT cont", "PLS NT mixed")
  return(results)
}

kfoldCV_reg3 <- function(dim_reduction=1,number_of_folds=5,binary_treatment=FALSE,ids){
  K <- number_of_folds
  a <- dim_reduction
  ids <- id
  # Create empty array to be filled by results:
  # 1st dimension: prediction methods (3 - linear model, non-parametric, inverse regression)
  # 2nd dimension: treatment and reduction methods (10)
  # 3rd dimension: number of folds
  Kresults <- array(NA, dim=c(K,10,2))
  for(i in 1:K){
    flag <- paste0("Running the ",i,"-th fold",sep="")
    print(flag)
    ytrain <- y[id!=i]
    ytest <- y[id==i]
    Xtrain <- X[id!=i,]
    Xtest <- X[id==i,]
    if(binary_treatment){
      #Projections with treatment of binary variables
      projections <- pca.pls(ytrain,ytest,Xtrain,Xtest,col.m,col.n,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }else{#Projections with NO treatment of binary variables
      projections <- pca.pls_binary(ytrain,ytest,Xtrain,Xtest,col.m,col.n,col.b,A=a)
      proj_train <- c(projections$pca_projections_train,projections$pls_projections_train)
      proj_test <- c(projections$pca_projections_test,projections$pls_projections_test)
    }
    for (pr in 1:length(proj_test)) {
      # Linear Regression
      fit <- linear_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MSE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 1] <- MSE_test
      # Non-Parametric Regression
      #fit <- np_model(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      #MSE_test <-  fit$MSE_TEST
      #Kresults[i , pr , 2] <- MSE_test
      # Inverse Regression
      fit <- inverse_reg(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]])
      # MCE
      MSE_test <-  fit$MSE_TEST
      Kresults[i , pr , 2] <- MSE_test
    }
  }
  results = colMeans(Kresults, dim=1)
  colnames(results) <- c("Linear Reg","Non-Parametric Reg")#, "Inverse Reg")
  rownames(results) <- c("PCA NM","PCA MNMC", "PCA OHE", "PCA NT cont", "PCA NT mixed",
                         "PLS NM","PLS MNMC", "PLS OHE", "PLS NT cont", "PLS NT mixed")
  return(results)
}
