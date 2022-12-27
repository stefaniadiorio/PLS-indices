p=50 # dim X
dx=3 # dimension de la regresion
q=15 # cantidad de no metricas
n=100

source("gen_ordinals.R")
source("data_processing/treatment.R")
source("reduction/pca.pls.R")
source("prediction_methods/prediction.R")

# Correlation Matrix
set.seed(2)
CorrMat = randcorr(p)
# Diagonal matrix of standard deviations

#D = diag(c(rep(1,(p-q)),runif(q, 0.2, 2)))
D=diag(rep(1,p))
# Covariance matrix
Sigma = D %*% CorrMat %*% D
mu = rep(0,p)

trueW=eigen(Sigma)$vectors[,c(7,8,13)]
#trueW=eigen(Sigma)$vectors[,c(2, 8 )]
a=3
M = 100 # number of iterations

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

for (s in 1:M){
  
  Z=mvrnorm(n, mu , Sigma)
  
  # outcome variable
  MM=  scale(Z,center=TRUE,scale=FALSE)%*%trueW
  
  
  y=rnorm(n, ((MM[,1]^(2)+sin(MM[,2])+MM[,3]^(2))) ,0.01)
  #y=(y+9)^(1/2)
  
  #y=rnorm(n, (MM[,1]+MM[,2])^(2) ,.1)
  ytrain=y[1:70]
  ytest=y[71:100]
  
  
  Xnm=Z[,1:q] # no metricas
  
  Xm=Z[,(q+1):(p)] # metricas
  
  
  # "traer" la función create_plist ######################################################### ¡!
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
  
  
  Xtrain=X[1:70,]
  Xtest=X[71:100,]
  
  col.n=c(1:q)
  col.m=c(q+1:(p-q))
  
  
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
  
  
  print(c("Finish iteration number: ", as.character(s)))
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

