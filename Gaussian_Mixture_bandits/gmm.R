# this is a script for testing our method on bandits with GMM






######################################################
###             Generate Test Problems             ###
######################################################
C = 5
piC = runif(C)
piC = piC/sum(piC)
alpha = 1:5
l = sample(x=alpha,size=1000,prob=piC,replace=TRUE)
l=sort(l)
# now, generate the arms according to a N(mu,1) distribution
mu = l + 0.25*rnorm(1000)
s <-list()
x<-vector()
# then, generate some samples from the each of the arm, then proceed to estimate the mixture model
for(i in 1:1000){
  s[[i]]<-vector()
  for(t in 1:1000){
    s[[i]][t]<-mu[i]+0.5*rnorm(1)
  }
  x<-c(x,s[[i]])
}
# from these data, we want to perform some estimation of the mixture model
avg<-vector()
for(i in 1:1000){
  avg[i]<-mean(s[[i]])
}

# write a function which could estimate both the arm parameters and the cluster parameters