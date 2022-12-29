### k-fold CV for classification ###
# Inputs:
# - dim_reduction (integer): number of directions PCA and PLS
# - number_of_folds (integer): number of K-folds
# - binary_treatment (logical): if binary variables are included or not in non-metric variables treatment
kfoldCV_clf <- function(dim_reduction=1,number_of_folds=5,binary_treatment=FALSE, ids){
  K <- number_of_folds
  a <- dim_reduction
  ids <- id
  # Create empty array to be filled by results:
  # 1st dimension: prediction methods (4 - log_regression, lda, qda, inverse regression)
  # 2nd dimension: treatment and reduction methods (10)
  # 3rd dimension: number of folds
  Kresults <- array(NA, dim=c(K,10,4))
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
      # Logistic Regression
      fit <- logistic(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]],c=0.5)
      # AUC
      preds <-  fit$Preds
      y_pred <-  preds$y_pred_test
      Kresults[i , pr , 1] <- AUC(y_pred,ytest)
      # LDA
      fit <- lda_class(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]],c=0.5)
      # AUC
      preds <-  fit$Preds
      y_pred <-  preds$y_pred_test
      Kresults[i , pr , 2] <- AUC(y_pred,ytest)
      # QDA
      fit <- qda_class(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]],c=0.5)
      # AUC
      preds <-  fit$Preds
      y_pred <-  preds$y_pred_test
      Kresults[i , pr , 3] <- AUC(y_pred,ytest)
      # Inverse Regression
      fit <- reg_inv(ytrain,ytest,Pr_train=proj_train[[pr]],Pr_test=proj_test[[pr]],c=0.5)
      # AUC
      preds <-  fit$Preds
      y_pred <-  preds$y_pred_test
      Kresults[i , pr , 4] <- AUC(y_pred,ytest)
    }
  }
  results = colMeans(Kresults, dim=1)
  colnames(results) <- c("Logistic Reg","LDA","QDA", "Inverse Reg")
  rownames(results) <- c("PCA NM","PCA MNMC", "PCA OHE", "PCA NT mixed", "PCA NT cont",
                         "PLS NM","PLS MNMC", "PLS OHE", "PLS NT mixed", "PLS NT cont")
  return(results)
}


