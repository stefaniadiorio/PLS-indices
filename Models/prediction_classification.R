### Prediction Methods ###
# INPUTS:
# - ytrain
# - ytest
# - Pr_train: projections over train set
# - Pr_test: projections over test set using loadings estimated with train set
# - c: cut

# OUTPUTS:
# MSE: mean squared error
#library(pROC)
# Logistic Regression
logistic <- function(ytrain,ytest,Pr_train,Pr_test,c){
  # Model Fitting
  logit.fit<-glm(ytrain~.,data=data.frame(ytrain,Pr_train),family=binomial(link = "logit"))
    datos = data.frame(ytest,Pr_test)
    names(datos)<-names(data.frame(ytrain,Pr_train))
  # Predict
  y_pred_train <- predict.glm(logit.fit, newdata = data.frame(ytrain,Pr_train), type="response")
  y_pred_test <- predict.glm(logit.fit, newdata = datos, type="response")
  # AUC
  #auc_train <- as.numeric(roc(ytrain~y_pred_train)$auc)
  #auc_test <- as.numeric(roc(ytest~y_pred_test)$auc)
  #AUC=list(auc_train=auc_train,auc_test=auc_test)
  # Mean Square Error
#  MSE_TRAIN <- mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
#  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
#  MSE=list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST)
  # Preds
  Preds=list(y_pred_train=y_pred_train, y_pred_test=y_pred_test)
  # Classify
#  y_pred_train[y_pred_train>=c]=1; y_pred_train[y_pred_train<c]=0
#  y_pred_test[y_pred_test>=c]=1; y_pred_test[y_pred_test<c]=0
  # MCE 
#  MCE_TRAIN=1-sum(y_pred_train==ytrain)/length(ytrain)
#  MCE_TEST=1-sum(y_pred_test==ytest)/length(ytest)
#  MCE=list(MCE_TRAIN=MCE_TRAIN, MCE_TEST=MCE_TEST)
  return(list(Preds=Preds))
}
# Linear Discriminant Analysis
lda_class <- function(ytrain,ytest,Pr_train,Pr_test,c){
  # Model Fitting
  lda.fit<-lda((ytrain)~.,data=data.frame(ytrain,Pr_train))
    datos = data.frame(ytest,Pr_test)
    names(datos)<-names(data.frame(ytrain,Pr_train))
  # Predict
  y_pred_train<-predict(lda.fit,newdata= data.frame(ytrain,Pr_train), type="response")$posterior[,2]
  y_pred_test<-predict(lda.fit,newdata= datos, type="response")$posterior[,2]
  # AUC
  #auc_train <- as.numeric(roc(ytrain~y_pred_train)$auc)
  #auc_test <- as.numeric(roc(ytest~y_pred_test)$auc)
  #AUC=list(auc_train=auc_train,auc_test=auc_test)
  # Mean Square Error
#  MSE_TRAIN <- mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
#  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
#  MSE=list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST)
  # Preds
  Preds=list(y_pred_train=y_pred_train, y_pred_test=y_pred_test)
  # Classify
#  y_pred_train[y_pred_train>=c]=1; y_pred_train[y_pred_train<c]=0
#  y_pred_test[y_pred_test>=c]=1; y_pred_test[y_pred_test<c]=0
  # MCE 
#  MCE_TRAIN=1-sum(y_pred_train==ytrain)/length(ytrain)
#  MCE_TEST=1-sum(y_pred_test==ytest)/length(ytest)
#  MCE=list(MCE_TRAIN=MCE_TRAIN, MCE_TEST=MCE_TEST)
  return(list(Preds=Preds))
}
# Quadratic Discriminant Analysis
qda_class <- function(ytrain,ytest,Pr_train,Pr_test,c){
  # Model Fitting
  qda.fit<-qda(ytrain~.,data=data.frame(ytrain,Pr_train))
  datos = data.frame(ytest,Pr_test)
  names(datos)<-names(data.frame(ytrain,Pr_train))
  # Predict
  y_pred_train<-predict(qda.fit,newdata= data.frame(ytrain,Pr_train), type="response")$posterior[,2]
  y_pred_test<-predict(qda.fit,newdata= datos, type="response")$posterior[,2]
  # AUC
  #auc_train <- as.numeric(roc(ytrain~y_pred_train)$auc)
  #auc_test <- as.numeric(roc(ytest~y_pred_test)$auc)
  #AUC=list(auc_train=auc_train,auc_test=auc_test)
  # Mean Square Error
#  MSE_TRAIN <- mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
#  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
#  MSE=list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST)
  # Preds
  Preds=list(y_pred_train=y_pred_train, y_pred_test=y_pred_test)
  # Classify
#  y_pred_train[y_pred_train>=c]=1; y_pred_train[y_pred_train<c]=0
#  y_pred_test[y_pred_test>=c]=1; y_pred_test[y_pred_test<c]=0
  # MCE 
#  MCE_TRAIN=1-sum(y_pred_train==ytrain)/length(ytrain)
#  MCE_TEST=1-sum(y_pred_test==ytest)/length(ytest)
#  MCE=list(MCE_TRAIN=MCE_TRAIN, MCE_TEST=MCE_TEST)
  return(list(Preds=Preds))
}

# Inverse Regression
reg_inv <- function(ytrain,ytest,Pr_train,Pr_test,c){
  # Model Fitting
  mod_inv=lm(Pr_train~ytrain)
  va=var(mod_inv$residuals)
  
  if (!is.matrix(mod_inv$fitted.values)){
    mod_inv$fitted.values <- matrix(mod_inv$fitted.values, ncol=1)
  }
  # Predict
  y_pred_train=inv_g(Pr_train,mod_inv$fitted.values,ytrain,va)
  y_pred_test=inv_g(Pr_test,mod_inv$fitted.values,ytrain,va)
  # AUC
  #auc_train <- as.numeric(roc(ytrain~y_pred_train)$auc)
  #auc_test <- as.numeric(roc(ytest~y_pred_test)$auc)
  #AUC=list(auc_train=auc_train,auc_test=auc_test)
  # Mean Square Error
#  MSE_TRAIN <- mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
#  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
#  MSE=list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST)
  # Preds
  Preds=list(y_pred_train=y_pred_train, y_pred_test=y_pred_test)
  # Classify
#  y_pred_train[y_pred_train>=c]=1; y_pred_train[y_pred_train<c]=0
#  y_pred_test[y_pred_test>=c]=1; y_pred_test[y_pred_test<c]=0
  # MCE 
#  MCE_TRAIN=1-sum(y_pred_train==ytrain)/length(ytrain)
#  MCE_TEST=1-sum(y_pred_test==ytest)/length(ytest)
#  MCE=list(MCE_TRAIN=MCE_TRAIN, MCE_TEST=MCE_TEST)
  return(list(Preds=Preds))
}


##### Auxiliar functions
inv_g<-function(z,prex,yy,va){
  w=NULL
  Rf=NULL
  for (h in 1:dim(z)[1]){
    for (i in 1:dim(prex)[1]){
      w[i]=exp(-1/2*t(z[h,]-prex[i,])%*%solve(va)%*%(z[h,]-prex[i,]))
    }
    R=sum(w*yy)/sum(w)
    Rf[h]=R
  }
  Rf
}

test_prediction <- function(clf_model, proj_train,proj_test,ytrain,ytest){
  if(clf_model == "Logistic.Reg"){
    fit <- logistic(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test,c=0.5)
    preds <-  fit$Preds
  }
  if(clf_model == "LDA"){
    fit <- lda_class(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test,c=0.5)
    preds <-  fit$Preds
  }
  if(clf_model == "QDA"){
    fit <- qda_class(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test,c=0.5)
    preds <-  fit$Preds

  }
  if(clf_model == "Inverse.Reg"){
    fit <- reg_inv(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test,c=0.5)
    preds <-  fit$Preds
  }
  return(preds)
}












