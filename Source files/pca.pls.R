source("data_processing/treatment.R")
source("reduction/pca-lili.R")
source("reduction/pls_lili.R")
pca.pls <- function(y_train,y_test,Xtrain,Xtest,col.m,col.n,A=1){
  datasets <- preprocessing(Xtrain,Xtest,col.m,col.n)
  ### Separate different treated datasets
  NMC <- datasets$NMC
  MNMC <- datasets$MNMC
  OHE <- datasets$OneHot
  NT <- datasets$NoTreat
  
  ### Train and Test data sets
    ## Normal Mean Coding
      Xnmc_train <- NMC[[1]]
      Xnmc_test <- NMC[[2]]
    ## Multivariate Normal Mean Coding
      Xmnmc_train <- MNMC[[1]]
      Xmnmc_test <- MNMC[[2]]
    ## One hot Encoding
      Xohe_train <- OHE[[1]]
      Xohe_test <- OHE[[2]]
    ## No-treatment
      Xnt_train <- NT[[1]]
      Xnt_test <- NT[[2]]
  
  
  ### Performs PCA and PLS on each treated Xtrain set and obtain loadings
  ### pca_reduction() and pls_reduction() returns loadings
  
  # 2.1.1 Normal mean coding (as continuous predictors)
    ## PCA
     w_nmc <- pca_lili(X=Xnmc_train, Xo=NA,A)
     w_nmc <- w_nmc$W_aut
     pls_nmc <- pls_lili(X=Xnmc_train,  Xo=NA,y=y_train, A)$W_aux
  
  # 2.1.2 Multivariate normal mean coding (as continuous predictors)
    ## PCA
      w_mnmc <- pca_lili(X=Xmnmc_train, Xo=NA, A)
      w_mnmc <- w_mnmc$W_aut
     pls_mnmc <- pls_lili(X=Xmnmc_train,  Xo=NA,y=y_train, A)$W_aux
  
  # 2.1.3 One hot encoding (as continuous predictors)
    ## PCA
      w_ohe_c <- pca_lili(X=Xohe_train, Xo=NA, A)
      w_ohe_c <- w_ohe_c$W_aut
     pls_ohe_c <- pls_lili(X=Xohe_train, Xo=NA, y=y_train, A)$W_aux
  
  ######## MIXTA continuas y ordinales...
      
      # 2.1.4 No treatment (mixed continuous and ordinals)
      ## PCA
      w_nt_m <- pca_lili(X=Xnt_train[,col.m], Xo=Xnt_train[,col.n], A)
      w_nt_m <- w_nt_m$W_aut
      pls_aux <- pls_lili(X=Xnt_train[,col.m], Xo=Xnt_train[,col.n], y=y_train, A)$W_aux
      
      npm=c(col.n,col.m)
      names(npm)<-seq(npm)
      npm=sort(npm)
      new_order=as.numeric(names(npm))
      pls_nt_m = pls_aux[new_order,]
      
  # 2.1.5 No treatment (as continuous predictors)
    ##PCA
      w_nt_c <- pca_lili(X=Xnt_train, Xo=NA,A)
      w_nt_c <- w_nt_c$W_aut
      pls_nt_c <- pls_lili(X=Xnt_train,  Xo=NA,y=y_train, A)$W_aux
  

  
  # Projections PCA
  pca_projections_train <- list(#Projections over Xtrain
              Pr.nmc_train = Xnmc_train%*%w_nmc,
              Pr.mnmc_train = Xmnmc_train%*%w_mnmc,
              Pr.ohe_c_train = Xohe_train%*%w_ohe_c,
              Pr.nt_c_train = Xnt_train%*%w_nt_c,
              Pr.nt_m_train = Xnt_train%*%w_nt_m)
              #Projections over Xtest using loadings estimated with Xtrain
  pca_projections_test <- list(Pr.nmc_test = Xnmc_test%*%w_nmc,
              Pr.mnmc_test = Xmnmc_test%*%w_mnmc,
              Pr.ohe_c_test = Xohe_test%*%w_ohe_c,
              Pr.nt_c_test = Xtest%*%w_nt_c,
              Pr.nt_m_test = Xtest%*%w_nt_m)

  # Projections PLS
  pls_projections_train <- list(#Projections over Xtrain
              Pr.nmc_train = Xnmc_train%*%pls_nmc,
              Pr.mnmc_train = Xmnmc_train%*%pls_mnmc,
              Pr.ohe_c_train = Xohe_train%*%pls_ohe_c,
              Pr.nt_c_train = Xnt_train%*%pls_nt_c,
              Pr.nt_m_train = Xnt_train%*%pls_nt_m)
              #Projections over Xtest using loadings estimated with Xtrain
  pls_projections_test <- list(Pr.nmc_test = Xnmc_test%*%pls_nmc,
              Pr.mnmc_test = Xmnmc_test%*%pls_mnmc,
              Pr.ohe_c_test = Xohe_test%*%pls_ohe_c,
              Pr.nt_c_test = Xtest%*%pls_nt_c,
              Pr.nt_m_test = Xtest%*%pls_nt_m)
  
  return(list(pca_projections_train=pca_projections_train,
              pca_projections_test=pca_projections_test,
              pls_projections_train=pls_projections_train,
              pls_projections_test=pls_projections_test))
}

pca.pls_binary <- function(y_train,y_test,Xtrain,Xtest,col.m,col.n,col.b,A=1){
  datasets <- preprocessing_binary(Xtrain,Xtest,col.m,col.n,col.b)
  ### Separate different treated datasets
  NMC <- datasets$NMC
  MNMC <- datasets$MNMC
  OHE <- datasets$OneHot
  NT <- datasets$NoTreat
  
  ### Train and Test data sets
  ## Normal Mean Coding
  Xnmc_train <- NMC[[1]]
  Xnmc_test <- NMC[[2]]
  ## Multivariate Normal Mean Coding
  Xmnmc_train <- MNMC[[1]]
  Xmnmc_test <- MNMC[[2]]
  ## One hot Encoding
  Xohe_train <- OHE[[1]]
  Xohe_test <- OHE[[2]]
  ## No-treatment
  Xnt_train <- NT[[1]]
  Xnt_test <- NT[[2]]
  
  
  ### Performs PCA and PLS on each treated Xtrain set and obtain loadings
  ### pca_reduction() and pls_reduction() returns loadings
  
  # 2.1.1 Normal mean coding (as continuous predictors)
  ## PCA
  w_nmc <- pca_lili(X=Xnmc_train, Xo=NA,A)
  w_nmc <- w_nmc$W_aut
  pls_nmc <- pls_lili(X=Xnmc_train,  Xo=NA,y=y_train, A)$W_aux
  
  # 2.1.2 Multivariate normal mean coding (as continuous predictors)
  ## PCA
  w_mnmc <- pca_lili(X=Xmnmc_train, Xo=NA, A)
  w_mnmc <- w_mnmc$W_aut
  pls_mnmc <- pls_lili(X=Xmnmc_train,  Xo=NA,y=y_train, A)$W_aux
  
  # 2.1.3 One hot encoding (as continuous predictors)
  ## PCA
  w_ohe_c <- pca_lili(X=Xohe_train, Xo=NA, A)
  w_ohe_c <- w_ohe_c$W_aut
  pls_ohe_c <- pls_lili(X=Xohe_train, Xo=NA, y=y_train, A)$W_aux
  
  ######## MIXTA continuas y ordinales...
  
  # 2.1.4 No treatment (mixed continuous and ordinals)
  ## PCA
  w_nt_m <- pca_lili(X=Xnt_train[,c(col.m,col.b)], Xo=Xnt_train[,col.n], A)
  w_nt_m <- w_nt_m$W_aut
  pls_aux <- pls_lili(X=Xnt_train[,c(col.m,col.b)], Xo=Xnt_train[,col.n], y=y_train, A)$W_aux
  
  npm=c(col.n,col.m,col.b)
  names(npm)<-seq(npm)
  npm=sort(npm)
  new_order=as.numeric(names(npm))
  pls_nt_m = pls_aux[new_order,]
  
  # 2.1.5 No treatment (as continuous predictors)
  ##PCA
  w_nt_c <- pca_lili(X=Xnt_train, Xo=NA,A)
  w_nt_c <- w_nt_c$W_aut
  pls_nt_c <- pls_lili(X=Xnt_train,  Xo=NA,y=y_train, A)$W_aux
  
  
  
  # Projections PCA
  pca_projections_train <- list(#Projections over Xtrain
    Pr.nmc_train = Xnmc_train%*%w_nmc,
    Pr.mnmc_train = Xmnmc_train%*%w_mnmc,
    Pr.ohe_c_train = Xohe_train%*%w_ohe_c,
    Pr.nt_c_train = Xnt_train%*%w_nt_c,
    Pr.nt_m_train = Xnt_train%*%w_nt_m)
  #Projections over Xtest using loadings estimated with Xtrain
  pca_projections_test <- list(Pr.nmc_test = Xnmc_test%*%w_nmc,
                               Pr.mnmc_test = Xmnmc_test%*%w_mnmc,
                               Pr.ohe_c_test = Xohe_test%*%w_ohe_c,
                               Pr.nt_c_test = Xtest%*%w_nt_c,
                               Pr.nt_m_test = Xtest%*%w_nt_m)
  
  # Projections PLS
  pls_projections_train <- list(#Projections over Xtrain
    Pr.nmc_train = Xnmc_train%*%pls_nmc,
    Pr.mnmc_train = Xmnmc_train%*%pls_mnmc,
    Pr.ohe_c_train = Xohe_train%*%pls_ohe_c,
    Pr.nt_c_train = Xnt_train%*%pls_nt_c,
    Pr.nt_m_train = Xnt_train%*%pls_nt_m)
  #Projections over Xtest using loadings estimated with Xtrain
  pls_projections_test <- list(Pr.nmc_test = Xnmc_test%*%pls_nmc,
                               Pr.mnmc_test = Xmnmc_test%*%pls_mnmc,
                               Pr.ohe_c_test = Xohe_test%*%pls_ohe_c,
                               Pr.nt_c_test = Xtest%*%pls_nt_c,
                               Pr.nt_m_test = Xtest%*%pls_nt_m)
  
  return(list(pca_projections_train=pca_projections_train,
              pca_projections_test=pca_projections_test,
              pls_projections_train=pls_projections_train,
              pls_projections_test=pls_projections_test))
} 

test_reduction <- function(reduction,Xtrain_t,Xtest_t,ytrain,A){
  if(binary_treatment == TRUE){
    if(reduction == "PCA"){
      if(treatment == "NT mixed"){
        w <- pca_lili(X=Xtrain_t[,col.m], Xo=Xtrain_t[,col.n], A)$W_aut
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }else{
        w <- pca_lili(X=Xtrain_t, Xo=NA,A)$W_aut
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }
    }
    if(reduction == "PLS"){
      if(treatment == "NT mixed"){
        w <- pls_lili(X=Xtrain_t[,col.m], Xo=Xtrain_t[,col.n], y=ytrain, A)$W_aux
        npm=c(col.n,col.m)
        names(npm)<-seq(npm)
        npm=sort(npm)
        new_order=as.numeric(names(npm))
        w = w[new_order,]
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }else{
        w <- pls_lili(X=Xtrain_t,  Xo=NA, y=ytrain, A)$W_aux
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }
    }
  }else{
    if(reduction == "PCA"){
      if(treatment == "NT mixed"){
        w <- pca_lili(X=Xtrain_t[,c(col.m,col.b)], Xo=Xtrain_t[,col.n], A)$W_aut
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }else{
        w <- pca_lili(X=Xtrain_t, Xo=NA,A)$W_aut
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }
    }
    if(reduction == "PLS"){
      if(treatment == "NT mixed"){
        w <- pls_lili(X=Xtrain_t[,c(col.m,col.b)], Xo=Xtrain_t[,col.n], y=ytrain, A)$W_aux
        npm=c(col.n,col.m,col.b)
        names(npm)<-seq(npm)
        npm=sort(npm)
        new_order=as.numeric(names(npm))
        w = w[new_order,]
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }else{
        w <- pls_lili(X=Xtrain_t,  Xo=NA,y=ytrain, A)$W_aux
        proj_train <- Xtrain_t%*%w
        proj_test <- Xtest_t%*%w
      }
    }
  }
  return(list(proj_train,proj_test))
}



