### Simulations ###

library(MLmetrics)
source("data_generation/data_generation_clf.R")
source("data_processing/treatment.R")
source("reduction/pca.pls.R")
source("prediction_methods/prediction_classification.R")

### Settings
n_train = 1000 #1000
n_test = 500 #500
k = 50 
pnm = 0.8 # proportion of non-metric variables
dgp = 4

M = 100 # number of repeticiones: usar 100

data_test = call_dgp(n_test,k,pnm,dgp)
col.m <- data_test$col.m
col.n <- data_test$col.n
Xtest <- data_test$X
#ytest <- data_test$y

#armar la variable binaria
ytest_clf = rep(0, nrow(data_test$y))
ytest_clf[data_test$y > median(data_test$y)] <- 1
ytest <- ytest_clf

a = 2

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

for (i in 1:M){
  data_train = call_dgp(n_train,k,pnm,dgp)
  #data_test = call_dgp(n_test,k,pnm,dgp)
  Xtrain = data_train$X
  #ytrain <- data_train$y
  
  #armar la variable binaria
  ytrain_clf = rep(0, nrow(data_train$y))
  ytrain_clf[data_train$y > median(data_train$y)] <- 1
  ytrain <- ytrain_clf
  
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
  
  
  print(c("Finish iteration number: ", as.character(i)))
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
