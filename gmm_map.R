# this is for the MAP estimation of gaussian mixture models


# we start by generating some initial configuration to feed into the EM algorithm
sigma <- vector()
for(i in 1:length(s)){
  sigma[i] <- sqrt(var(s[[i]]))
}
# even with n_k = 1000, the estimation of the variance can still be a little off

# first, cluster all the samples, and use this as an initial configuration for MAP
p = rep(0.2,times=5)
v = rep(var(avg)/5,times = 5)
mu = rep(min(avg),times=5) + (max(avg)-min(avg))/6*(1:5)
w<-matrix(0,nrow = 5, ncol=1000)
for(k in 1:1000){
  for(j in 1:5){
    w[j,k] = p[j]*dnorm(avg[k],mean = mu[j], sd=v[j])
  }
  w[,k] = w[,k]/sum(w[,k])
}
# now, we can proceed to the iterative updates
# how to determine convergence? Let's see
d = 1
while(d>10^{-6}){
  # update the cluster parameters
  p0 = p
  for(j in 1:5){
    p[j] = sum(w[j,])/1000
    mu[j] = sum(w[j,]*avg)/sum(w[j,])
    v[j] = sqrt( sum( w[j,]*(avg - mu[j])^2 )/sum(w[j,]))
  } 
  # update the weights
  for(k in 1:1000){
    for(j in 1:5){
      w[j,k] = p[j]*dnorm(avg[k],mean=mu[j],sd=v[j])
    }
    w[,k] = w[,k]/sum(w[,k])
  }
  d = sum((p0-p)^2)
}


###################################
#    the estimation function      #
###################################

# initilization function, generate clusters from only the average of the samples
K = 1000 # number of arms
C = 5 # number of prior clusters
nk = rep(1000,times =1000 )

problem_set <-list(C,K,nk)
reward_sample <- s
alpha = 1:C
beta = rep(0.25,times = C)
p = rep(1/C,times=C)
para_cluster<-list(alpha,beta,p); names(para_cluster) = c('alpha','beta','p')
mu = avg
sigma = rep(0.5,times = K) 
para_arm = list(mu,sigma); names(para_arm) = c('mu','sigma')
para<-list(para_cluster,para_arm)
# with these data and initial set up, we move onto the optimization scheme

#
# this is the function we can use for the estimation problem
#

gmm_opt<-function(reward_sample,para,problem_set){
  
  w<-matrix(0,nrow=problem_set[[1]],ncol=problem_set[[2]])
  d = 1
  alpha = para[[1]][[1]]; beta = para[[1]][[2]]; p = para[[1]][[3]]; mu = para[[2]][[1]]; sigma = para[[2]][[2]];
  while(d>10^(-10)){
    # keeping track of convergence steps
    mu0 = mu
    # first, update the weights
      for(k in 1:problem_set[[2]]){
        for(j in 1:problem_set[[1]]){
          w[j,k] = p[j]*dnorm(mu[k],mean=alpha[j],sd=beta[j])
        }
        w[,k] = w[,k]/sum(w[,k])
      }   
    # then, use the weights to update all cluster parameters
    for(j in 1:problem_set[[1]]){
      p[j] = sum(w[j,])/problem_set[[2]]
      alpha[j] = sum( w[j,]*mu )/(p[j]*problem_set[[2]])
      beta[j] = sqrt(  sum(w[j,]*(mu-alpha[j])^2)/(p[j]*problem_set[[2]])   )
    }
    # update all the arm parameters
    for(k in 1:problem_set[[2]]){
      sigma[k] = 1/problem_set[[3]][k]*(sum( (reward_sample[[k]] - mu[k] )^2  ))
      mu[k] = ( sum(w[,k]*alpha/beta^2) + 1/sigma[k]^2*sum( reward_sample[[k]] ) )/( sum( w[,k]/beta^2 ) + problem_set[[3]][k]/sigma[k]^2  )
    }
    d = sum( (mu-mu0)^2 )
  }
  para_result<-list( list(alpha,beta,p), list(mu,sigma)  )
  return(para_result)
}












