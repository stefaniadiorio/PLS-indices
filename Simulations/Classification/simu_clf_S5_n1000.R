p=50 # dim X
dx=3 # dimension de la regresion
q=15 # cantidad de no metricas
n=1000

library(MLmetrics)
source("gen_ordinals.R")
source("data_processing/treatment.R")
source("reduction/pca.pls.R")
source("prediction_methods/prediction_classification.R")


# Correlation Matrix
set.seed(2)
CorrMat = randcorr(p)
# Diagonal matrix of standard deviations

D=diag(rep(1,p))
# Covariance matrix
Sigma = D %*% CorrMat %*% D
mu = rep(0,p)

trueW=eigen(Sigma)$vectors[,c(7,8,13)]

a=3
M = 100 # number of iterations

AUC_train_pca_log <- c()
AUC_test_pca_log <- c()
AUC_train_pls_log <- c()
AUC_test_pls_log <- c()

AUC_train_pca_lda <- c()
AUC_test_pca_lda <- c()
AUC_train_pls_lda <- c()
AUC_test_pls_lda <- c()

AUC_train_pca_qda <- c()
AUC_test_pca_qda <- c()
AUC_train_pls_qda <- c()
AUC_test_pls_qda <- c()

for (s in 1:M){
  
  Z=mvrnorm(n, mu , Sigma)
  
  # outcome variable
  MM=scale(Z,center=TRUE,scale=FALSE)%*%trueW
  
  y=rnorm(n,MM[,1] +MM[,2]+MM[,3],.5)
  y=(y+9)^(1/2)
  ytrain=as.matrix(y[1:700])
  ytrain_clf = rep(0, nrow(ytrain))
  ytrain_clf[ytrain > median(ytrain)] <- 1
  ytrain <- ytrain_clf
  
  ytest=as.matrix(y[701:1000])
  ytest_clf = rep(0, nrow(ytest))
  ytest_clf[ytest > median(ytest)] <- 1
  ytest <- ytest_clf
  
  Xnm=Z[,1:q] # ordinales
  
  Xm=Z[,(q+1):(p)] # metricas
  
  # para ordinales ################ 
  plist <- create_plist(p=q,max_num_categories = 7, min_num_categories=5)
  
  # generalizar a p variables
  Xnm = cbind()
  for(i in 1:q){
    breaks_i = as.matrix(qnorm(plist[[i]], mean=mean(Z[,i]), sd=sd(Z[,i])))
    breaks_i=c(-Inf, breaks_i, Inf)
    Z_i = cut(Z[,i], breaks = breaks_i, include.lowest=TRUE, right=TRUE, labels=FALSE)
    Xnm = cbind(Xnm,Z_i)
  }
  
  X = cbind(Xnm,Xm)
  
  
  Xtrain=X[1:700,]
  Xtest=X[701:1000,]
  
  col.n=c(1:q)
  col.m=c(q+1:(p-q))
  
  
  projections <- pca.pls(ytrain,ytest,Xtrain,Xtest,col.m,col.n,A=a)
  pca_proj_train <- projections$pca_projections_train
  pca_proj_test <- projections$pca_projections_test
  pls_proj_train <- projections$pls_projections_train
  pls_proj_test <- projections$pls_projections_test
  
  auc_train_pca_log <- c()
  auc_test_pca_log <- c()
  auc_train_pls_log <- c()
  auc_test_pls_log <- c()
  
  auc_train_pca_lda <- c()
  auc_test_pca_lda <- c()
  auc_train_pls_lda <- c()
  auc_test_pls_lda <- c()
  
  auc_train_pca_qda <- c()
  auc_test_pca_qda <- c()
  auc_train_pls_qda <- c()
  auc_test_pls_qda <- c()
  
  ### Prediction
  # Logistic regression (log), Linear Discriminant Analysis (lda), 
  #Quadratic Discriminant Analysis (qda)
  for (pr in 1:length(pca_proj_train)) {
    
    # Logistic regression
    #PCA
    Preds <- logistic(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]], c=0.5)
    auc_train_pca_log <- c(auc_train_pca_log, AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pca_log <- c(auc_test_pca_log, AUC(Preds[[1]]$y_pred_test, ytest))
    #PLS
    Preds <- logistic(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]], c=0.5)
    auc_train_pls_log <- c(auc_train_pls_log, AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pls_log <- c(auc_test_pls_log, AUC(Preds[[1]]$y_pred_test, ytest))
    
    # Linear Discriminant Analysis
    #PCA
    Preds<- lda_class(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]], c=0.5)
    auc_train_pca_lda <- c(auc_train_pca_lda, AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pca_lda <- c(auc_test_pca_lda, AUC(Preds[[1]]$y_pred_test, ytest))
    #PLS
    Preds <- lda_class(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]], c=0.5)
    auc_train_pls_lda <- c(auc_train_pls_lda, AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pls_lda <- c(auc_test_pls_lda, AUC(Preds[[1]]$y_pred_test, ytest))
    
    # Quadratic Discriminant Analysis
    #PCA
    Preds <- qda_class(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]], c=0.5)
    auc_train_pca_qda <- c(auc_train_pca_qda, AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pca_qda <- c(auc_test_pca_qda, AUC(Preds[[1]]$y_pred_test, ytest))
    #PLS
    Preds <- qda_class(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]], c=0.5)
    auc_train_pls_qda <- c(auc_train_pls_qda,AUC(Preds[[1]]$y_pred_train, ytrain))
    auc_test_pls_qda <- c(auc_test_pls_qda,AUC(Preds[[1]]$y_pred_test, ytest))
  }
  
  AUC_train_pca_log <- rbind(AUC_train_pca_log,auc_train_pca_log)
  AUC_test_pca_log <- rbind(AUC_test_pca_log,auc_test_pca_log)
  AUC_train_pls_log <- rbind(AUC_train_pls_log,auc_train_pls_log)
  AUC_test_pls_log <- rbind(AUC_test_pls_log,auc_test_pls_log)
  
  AUC_train_pca_lda <- rbind(AUC_train_pca_lda,auc_train_pca_lda)
  AUC_test_pca_lda <- rbind(AUC_test_pca_lda,auc_test_pca_lda)
  AUC_train_pls_lda <- rbind(AUC_train_pls_lda,auc_train_pls_lda)
  AUC_test_pls_lda <- rbind(AUC_test_pls_lda,auc_test_pls_lda)
  
  AUC_train_pca_qda <- rbind(AUC_train_pca_qda,auc_train_pca_qda)
  AUC_test_pca_qda <- rbind(AUC_test_pca_qda,auc_test_pca_qda)
  AUC_train_pls_qda <- rbind(AUC_train_pls_qda,auc_train_pls_qda)
  AUC_test_pls_qda <- rbind(AUC_test_pls_qda,auc_test_pls_qda)
  
  
  print(c("Finish iteration number: ", as.character(s)))
}


# Dataframe Train
AUC_train <- rbind(AUC_train_pca_log,AUC_train_pls_log,AUC_train_pca_lda,AUC_train_pls_lda,AUC_train_pca_qda,AUC_train_pls_qda)
colnames(AUC_train) <- c("NMC","MNMC","OHE","NT_m","NT_c")
Reduction <-c()
for (string in c(row.names(AUC_train))){splits <- strsplit(string,"_");Reduction <- c(Reduction,splits[[1]][3])}
AUC_train <- cbind(AUC_train,Reduction)
Prediction_Method <- c()
for (string in c(row.names(AUC_train))){splits <- strsplit(string,"_");Prediction_Method <- c(Prediction_Method,splits[[1]][4])}
AUC_train <- cbind(AUC_train,Prediction_Method)
row.names(AUC_train) <- NULL
AUC_train <- as.data.frame(AUC_train)
AUC_train[,1:5] <- lapply(AUC_train[,1:5], function(x) as.numeric(as.character(x)))

# Dataframe Test 
AUC_test <- rbind(AUC_test_pca_log,AUC_test_pls_log,AUC_test_pca_lda,AUC_test_pls_lda,AUC_test_pca_qda,AUC_test_pls_qda)
colnames(AUC_test) <- c("NMC","MNMC","OHE","NT_m","NT_c")
Reduction <-c()
for (string in c(row.names(AUC_test))){splits <- strsplit(string,"_");Reduction <- c(Reduction,splits[[1]][3])}
AUC_test <- cbind(AUC_test,Reduction)
Prediction_Method <- c()
for (string in c(row.names(AUC_test))){splits <- strsplit(string,"_");Prediction_Method <- c(Prediction_Method,splits[[1]][4])}
AUC_test <- cbind(AUC_test,Prediction_Method)
row.names(AUC_test) <- NULL
AUC_test <- as.data.frame(AUC_test)
AUC_test[,1:5] <- lapply(AUC_test[,1:5], function(x) as.numeric(as.character(x)))

