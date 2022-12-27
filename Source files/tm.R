# tm

# tm returns the number of categories
# as a function of expected number of categories
# and number of non-metric variables.

# inputs
# enc: expected number of categories. It should be 3 or 7
# nnmv: number of non-metric variables

# outputs
# results: a vector containing number of categories

# write.dta(data.frame(mis(n=rpo=50, lambda=1)+2), "m3.dta")
# write.dta(data.frame(m=rpois(n=50, lambda=5)+2), "m7.dta")

tm<-function(enc, nnmv){
  if(enc==3){
    results<-data.matrix(read.dta("data_generation/m3.dta"))
  } else {
    results<-data.matrix(read.dta("data_generation/m7.dta"))
  }
  
  results<-as.numeric(results[1:nnmv])
  
  return(results)
}

tm_binary<-function(enc_bin, nbin){
  if(enc_bin==2){
    results<-data.matrix(rep(2,50))
  } else {
    results<-data.matrix(read.dta("data_generation/m7.dta"))
  }
  
  results<-as.numeric(results[1:nbin])
  
  return(results)
}


###
# tf.1 and tf.2

# These programs help me to input parameters as a function of dist.xi
# It will be used in efgs.

# inputs
# i: dist.xi value
# i==1->normal; i!=1->lognormal

# outputs
# results: mean and sigma^2 parameters of normal or lognormal

tf.1<-function(i){
  if(i==1){return(c(0,1))}
  if(i!=1){return(c(-1.435311, 1.554052))}
}

# inputs
# i: dist.xi value
# i==1->normal; i!=1->lognormal

# outputs
# results: mean and sigma^2 parameters of normal or lognormal

tf.2<-function(i){
  if(i==1){return(c(0,5))}
  if(i!=1){return(c(-0.6305921,  1.5540524))}
}

###

