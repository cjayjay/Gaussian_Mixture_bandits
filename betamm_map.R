# this script is for solving the beta mixture model problem

z = rbeta(10,shape1 = 5, shape2 = 99)
library(RPMM)

loglik_beta<-function(mu,x){
  return( sum( -dbeta(x,mu[1],mu[2],log=TRUE)  ) )
}

out<-optim(par=c(1,1), fn =loglik_beta, x = z, method="L-BFGS-B",lower=c(0,0) )

# the method of moment estimator 
beta_MOM<-function(z){
  t1 = mean(z)
  t2 = var(z)
  a = t1^2*(1-t1)/t2-t1
  b = t1*(1-t1)^2/t2+t1-1
  return(c(a,b))
}
# the MLE estimator
beta_MLE<-function(z,x0){
  loglik_beta<-function(mu,x){
    return(sum(-dbeta(x,mu[1],mu[2],log=TRUE)))
  }
  out<-optim(par=x0,fn=loglik_beta, x=z, method = "L-BFGS-B",lower=c(0,0) )
  return(out[[1]])
}
# the difference seems to be pretty small



##################################################
#    the mixture model optimization function     #
##################################################
