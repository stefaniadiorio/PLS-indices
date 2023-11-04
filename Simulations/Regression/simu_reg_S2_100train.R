### Simulations ###


source("data_generation/data_generation.R")
source("data_processing/treatment.R")
source("reduction/pca.pls.R")
source("prediction_methods/prediction.R")

### Settings
n_train = 100
n_test = 500
k = 50 
pnm = 0.8 # proportion of non-metric variables
dgp = 2 #in dgp-lili choose the weights equal 1

M = 100 # number of repeticiones

data_test = call_dgp(n_test,k,pnm,dgp)
col.m <- data_test$col.m
col.n <- data_test$col.n
Xtest <- data_test$X
ytest <- data_test$y
a = 1

MSE_train_pca_lm <- c()
MSE_test_pca_lm <- c()
MSE_train_pls_lm <- c()
MSE_test_pls_lm <- c()

MSE_train_pca_np <- c()
MSE_test_pca_np <- c()
MSE_train_pls_np <- c()
MSE_test_pls_np <- c()

MSE_train_pca_ir <- c()
MSE_test_pca_ir <- c()
MSE_train_pls_ir <- c()
MSE_test_pls_ir <- c()

for (i in 1:M){
  data_train = call_dgp(n_train,k,pnm,dgp)
  data_test = call_dgp(n_test,k,pnm,dgp)
  Xtrain = data_train$X
  ytrain <- data_train$y
  
  projections <- pca.pls(ytrain,ytest,Xtrain,Xtest,col.m,col.n,A=a)
  pca_proj_train <- projections$pca_projections_train
  pca_proj_test <- projections$pca_projections_test
  pls_proj_train <- projections$pls_projections_train
  pls_proj_test <- projections$pls_projections_test
  
  mse_train_pca_lm <- c()
  mse_test_pca_lm <- c()
  mse_train_pls_lm <- c()
  mse_test_pls_lm <- c()
  
  mse_train_pca_np <- c()
  mse_test_pca_np <- c()
  mse_train_pls_np <- c()
  mse_test_pls_np <- c()
  
  mse_train_pca_ir <- c()
  mse_test_pca_ir <- c()
  mse_train_pls_ir <- c()
  mse_test_pls_ir <- c()
  
  ### Prediction
  # Linear model (lm), Non-Parametric Regression (np), Inverse Regression (ir)
  for (pr in 1:length(pca_proj_train)) {
    
    # Linear regression
    #PCA
    MSE <- linear_model(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]])
    mse_train_pca_lm <- c(mse_train_pca_lm,MSE$MSE_TRAIN)
    mse_test_pca_lm <- c(mse_test_pca_lm,MSE$MSE_TEST)
    #PLS
    MSE <- linear_model(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]])
    mse_train_pls_lm <- c(mse_train_pls_lm,MSE$MSE_TRAIN)
    mse_test_pls_lm <- c(mse_test_pls_lm,MSE$MSE_TEST)
    # Non-Parametric
    #PCA
    MSE <- np_model(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]])
    mse_train_pca_np <- c(mse_train_pca_np,MSE$MSE_TRAIN)
    mse_test_pca_np <- c(mse_test_pca_np,MSE$MSE_TEST)
    #PLS
    MSE <- np_model(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]])
    mse_train_pls_np <- c(mse_train_pls_np,MSE$MSE_TRAIN)
    mse_test_pls_np <- c(mse_test_pls_np,MSE$MSE_TEST)
    
    # Inverse Regression
    #PCA
    MSE <- inverse_reg(ytrain,ytest,Pr_train=pca_proj_train[[pr]],pca_proj_test[[pr]])
    mse_train_pca_ir <- c(mse_train_pca_ir,MSE$MSE_TRAIN)
    mse_test_pca_ir <- c(mse_test_pca_ir,MSE$MSE_TEST)
    #PLS
    MSE <- inverse_reg(ytrain,ytest,Pr_train=pls_proj_train[[pr]],pls_proj_test[[pr]])
    mse_train_pls_ir <- c(mse_train_pls_ir,MSE$MSE_TRAIN)
    mse_test_pls_ir <- c(mse_test_pls_ir,MSE$MSE_TEST)
  }
  
  MSE_train_pca_lm <- rbind(MSE_train_pca_lm,mse_train_pca_lm)
  MSE_test_pca_lm <- rbind(MSE_test_pca_lm,mse_test_pca_lm)
  MSE_train_pls_lm <- rbind(MSE_train_pls_lm,mse_train_pls_lm)
  MSE_test_pls_lm <- rbind(MSE_test_pls_lm,mse_test_pls_lm)
  
  MSE_train_pca_np <- rbind(MSE_train_pca_np,mse_train_pca_np)
  MSE_test_pca_np <- rbind(MSE_test_pca_np,mse_test_pca_np)
  MSE_train_pls_np <- rbind(MSE_train_pls_np,mse_train_pls_np)
  MSE_test_pls_np <- rbind(MSE_test_pls_np,mse_test_pls_np)
  
  MSE_train_pca_ir <- rbind(MSE_train_pca_ir,mse_train_pca_ir)
  MSE_test_pca_ir <- rbind(MSE_test_pca_ir,mse_test_pca_ir)
  MSE_train_pls_ir <- rbind(MSE_train_pls_ir,mse_train_pls_ir)
  MSE_test_pls_ir <- rbind(MSE_test_pls_ir,mse_test_pls_ir)
  
  
  print(c("Finish iteration number: ", as.character(i)))
}

# Dataframe Train
MSE_train <- rbind(MSE_train_pca_lm,MSE_train_pls_lm,MSE_train_pca_np,MSE_train_pls_np,MSE_train_pca_ir,MSE_train_pls_ir)
colnames(MSE_train) <- c("NMC","MNMC","OHE","NT_m","NT_c")
Reduction <-c()
for (string in c(row.names(MSE_train))){splits <- strsplit(string,"_");Reduction <- c(Reduction,splits[[1]][3])}
MSE_train <- cbind(MSE_train,Reduction)
Prediction_Method <- c()
for (string in c(row.names(MSE_train))){splits <- strsplit(string,"_");Prediction_Method <- c(Prediction_Method,splits[[1]][4])}
MSE_train <- cbind(MSE_train,Prediction_Method)
row.names(MSE_train) <- NULL
MSE_train <- as.data.frame(MSE_train)
MSE_train[,1:5] <- lapply(MSE_train[,1:5], function(x) as.numeric(as.character(x)))

# Dataframe Test 
MSE_test <- rbind(MSE_test_pca_lm,MSE_test_pls_lm,MSE_test_pca_np,MSE_test_pls_np,MSE_test_pca_ir,MSE_test_pls_ir)
colnames(MSE_test) <- c("NMC","MNMC","OHE","NT_m","NT_c")
Reduction <-c()
for (string in c(row.names(MSE_test))){splits <- strsplit(string,"_");Reduction <- c(Reduction,splits[[1]][3])}
MSE_test <- cbind(MSE_test,Reduction)
Prediction_Method <- c()
for (string in c(row.names(MSE_test))){splits <- strsplit(string,"_");Prediction_Method <- c(Prediction_Method,splits[[1]][4])}
MSE_test <- cbind(MSE_test,Prediction_Method)
row.names(MSE_test) <- NULL
MSE_test <- as.data.frame(MSE_test)
MSE_test[,1:5] <- lapply(MSE_test[,1:5], function(x) as.numeric(as.character(x)))



