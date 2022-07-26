rm(list=ls())

dev.off()

par(mfrow=c(1,1))

setwd("H:/Post doc ENSAI/PNAS")

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

####################################

#n=1000

#Z=sample(x = 0:100, n, replace = T, prob = ddgpd(0:100,0.5,1))

####################################

Z=c(580,68,39,275,46,61,40,38,73,64,134,21,1,1,38,50,9,26,30,197,57,12,22,113,26,2,7,19,68,39,46,33,87,85,39,24,19,49,55,31,48,132,63,27,17,54,46,12,13,22,9,44,36,49,47,61,17,12,23,12,4,13,13,14,13,4,16,6,6,6,17,11,66,15,4,3,3,2,2,2,13,18,33,17,31,20,18,15,38,7,13,21,3,18,22,24,84,39,22,20,18,67,66,27,18,34,41,17,20,29,35,17,12)

ntot=length(Z)

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
  
  seuil=which(vecu>=22)[1]
  
  xihat[j]=vecgamma[seuil]
  
  print(j*100/nrep)
  
  j=j+1
}

#save.image(file="robust_kerala_13_07.Rdata")