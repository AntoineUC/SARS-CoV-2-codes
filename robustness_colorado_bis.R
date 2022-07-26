rm(list=ls())

library(maps)

library(evt0)

library(Expectrem)

dev.off()

#par(mfrow=c(1,1))

setwd("H:/Post doc ENSAI/PNAS")

Z=read.table(file="colorado_new.txt")

Z=Z[,1]

Zfull=Z

ntot=length(Zfull)

nrep=10000

xihat=rep(0,nrep)

j=1

while(j<=nrep){
  
  #Ycalifornia=trunc(Zfull*runif(length(Zfull),0.9,1.1))+1
  
  Ycalifornia=trunc(Zfull*(1+rnorm(length(Zfull),sd=0.05)))+1
  
  n=length(Ycalifornia)
  
  #################discrete GPD###############
  
  Z=Ycalifornia
  
  ddgpd=function(k,igamma,sigma){
    return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
  }
  
  n=length(Z)
  
  Ztot=Z
  
  ###################################
  
  llfunction=function(theta){
    igamma=theta[1]
    sigma=theta[2]
    return(-sum(log((ddgpd(Z,igamma,sigma)))))
  }
  
  llfunction2=function(igamma,sigma){
    return(-sum(log(((1+igamma*Z/sigma)^(-1/igamma)-(1+igamma*(Z+1)/sigma)^(-1/igamma)))))
  }
  
  k=1
  
  #vecu=1:20
  
  vecu=unique(sort(Ztot))
  
  vecgamma=rep(0,length(vecu))
  
  vecgammainf=rep(0,length(vecu))
  
  vecgammasup=rep(0,length(vecu))
  
  vecsigma=rep(0,length(vecu))
  
  vecmoment=rep(0,length(vecu))
  
  vecsigmainf=rep(0,length(vecu))
  
  vecsigmasup=rep(0,length(vecu))
  
  while(k<=length(vecu)){
    
    Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
    
    nk=length(Z)
    
    thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par
    
    vecgamma[k]=thetahat[1]
    
    vecsigma[k]=thetahat[2]
    
    k=k+1
    
  }
  
  seuil=which(vecu>=27)[1]
  
  xihat[j]=vecgamma[seuil]
  
  print(j*100/nrep)
  
  j=j+1
}

#save.image(file="robust_colorado_13_07.Rdata")