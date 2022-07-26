rm(list=ls())
set.seed(1342)

# One single replication

n=10000
U=runif(n)

xi=1/2
rho=-1

X=(U^(rho)-1)^(-xi/rho)
Xround=trunc(X)+1

eps=10^(-1)
thresh=seq(eps,max(X)+1,by=eps)

Hill=rep(NA,length(thresh))
Hill2=rep(NA,length(thresh))

for (i in (1:length(thresh))){
  Hill[i]=mean(log(X[X>thresh[i]]/thresh[i]))
  Hill2[i]=mean(log(Xround[Xround>thresh[i]]/thresh[i]))
}

par(mfrow=c(1,2))
plot(thresh,Hill,type="l",xlim=c(0,10),ylim=c(0,1),xlab="Threshold",ylab="Estimate",main="Hill plot - Burr data")
abline(h=xi,lty=3)
plot(thresh,Hill2,type="l",xlim=c(0,10),ylim=c(0,1),xlab="Threshold",ylab="Estimate",main="Hill plot - Rounded Burr data")
abline(h=xi,lty=3)

# Repeating the experiment

N=1000
xi=1/2
rho=-1

eps=10^(-1)
thresh=seq(eps,10,by=eps)

Hillav=rep(0,length(thresh))
Hillav2=rep(0,length(thresh))

n=10000

for (j in 1:N){
  
  U=runif(n) # Uniform
  X=(U^(rho)-1)^(-xi/rho) # Burr
  Xround=trunc(X)+1 # Rounded Burr
  
  Hill=rep(NA,length(thresh))
  Hill2=rep(NA,length(thresh))
  
  for (i in (1:length(thresh))){
    Hill[i]=mean(log(X[X>thresh[i]]/thresh[i])) # Hill original data
    Hill2[i]=mean(log(Xround[Xround>thresh[i]]/thresh[i])) # Hill rounded data
    Hillav[i]=Hillav[i]+Hill[i]/N # Calculating the average Hill estimate over all replications
    Hillav2[i]=Hillav2[i]+Hill2[i]/N # Calculating the average Hill estimate over all replications (rounded data)
  }
  
}

par(mfrow=c(1,2))
plot(thresh,Hillav,type="l",xlim=c(0,10),ylim=c(0,1),xlab="Threshold",ylab="Average estimate",main="Average Hill plot")
abline(h=xi,lty=2)
plot(thresh,Hillav2,type="l",xlim=c(0,10),ylim=c(0,1),xlab="Threshold",ylab="Average estimate",main="Average Hill plot - Rounded data")
abline(h=xi,lty=2)