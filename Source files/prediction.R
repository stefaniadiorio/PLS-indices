### Prediction Methods ###
# INPUTS:
# - ytrain
# - ytest
# - Pr_train: projections over train set
# - Pr_test: projections over test set using loadings estimated with train set

# OUTPUTS:
# MSE: mean squared error
require(np)
# Linear Model
linear_model <- function(ytrain,ytest,Pr_train,Pr_test){
  lm.fit <- lm("ytrain~.",data=data.frame(ytrain,Pr_train), na.action="na.exclude")
  y_pred_train <- predict.lm(lm.fit, newdata = data.frame(ytrain,Pr_train))
  datos=data.frame(ytest,Pr_test)
  names(datos) <- names(data.frame(ytrain,Pr_train))
  y_pred_test <- predict.lm(lm.fit, newdata = datos)
  MSE_TRAIN<-mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
  return(list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST))
}

# Non-parametric
np_model <- function(ytrain,ytest,Pr_train,Pr_test){
 
  
  #bwnp_prueba<-np::npregbw(xdat=Pr_train,ydat=as.vector(ytrain))
  
  datas = data.frame(y = ytrain,X=Pr_train)
  A=ncol(Pr_train)
  model = "y~"
  if (A==1){
    model = paste(model,"X",sep="")
  }else{for (i in 1:A){
    if (i < A){
      x = paste("X.",i,sep = "")
      model = paste(model,x,"+",sep="")
    }else{
      x = paste("X.",i,sep = "")
      model = paste(model,x,sep="")
    }
  }}
  
 
  bwnp_prueba<-np::npregbw(as.formula(model), data=datas)

  
  bw <- np::npreg(bwnp_prueba, na.action="na.exclude")
  #bw <- npreg(tydat=as.vector(y), txdat=datos[,-1], regtype="ll", bwmethod="cv.aic") 
  #plot(y,predict(bw))
  datatrain<-data.frame(X=Pr_train)
  datatest=data.frame(X=Pr_test)
  y_pred_train <- predict(bw, newdata=datatrain, type="response")
  y_pred_test <- predict(bw, newdata=datatest, type="response")
   
  
  
 
 
  MSE_TRAIN <- mean(c(ytrain-y_pred_train)^2, na.rm=TRUE)
  MSE_TEST <- mean(c(ytest-y_pred_test)^2, na.rm=TRUE)
  return(list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST))
}

# Inverse regression
inverse_reg <- function(ytrain,ytest,Pr_train,Pr_test){
  
  aux=(abs(ytrain))^(1/2)
 
  
  mod_inv=lm(Pr_train~ytrain+aux+I(ytrain^2)+I(ytrain^3)+I(abs(ytrain)^(1/3)), na.action="na.exclude")
  va=var(mod_inv$residuals)
  
  if (!is.matrix(mod_inv$fitted.values)){
    mod_inv$fitted.values <- matrix(mod_inv$fitted.values, ncol=1)
  }
  
  pred_train=inv_g(Pr_train,mod_inv$fitted.values,ytrain,va)
  pred_test=inv_g(Pr_test,mod_inv$fitted.values,ytrain,va) 
  MSE_TRAIN <- mean(c(ytrain-pred_train)^2, na.rm=TRUE)
  MSE_TEST <- mean(c(ytest-pred_test)^2, na.rm=TRUE)
  return(list(MSE_TRAIN=MSE_TRAIN, MSE_TEST=MSE_TEST))
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



test_prediction <- function(reg_model, proj_train,proj_test,ytrain,ytest){
  if(reg_model == "Linear.Reg"){
    fit <- linear_model(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test)
    MSE_test <-  fit$MSE_TEST
  }
  if(reg_model == "Non.Parametric.Reg"){
    fit <- np_model(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test)
    MSE_test <-  fit$MSE_TEST
  }
  if(reg_model == "Inverse.Reg"){
    fit <- inverse_reg(ytrain,ytest,Pr_train=proj_train,Pr_test=proj_test)
    MSE_test <-  fit$MSE_TEST
  }
  return(MSE_test)
}


