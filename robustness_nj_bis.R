rm(list=ls())

library(maps)

library(evt0)

library(Expectrem)

dev.off()

#par(mfrow=c(1,1))

setwd("H:/Post doc ENSAI/PNAS")

country=read.table(file="country cluster.txt",sep="\n")

country=country[,1]

sum(country=="United States")

usindex=which(country=="United States")

total.cases=read.table(file="total cases.txt",header=F,sep="\n")

total.cases=total.cases[,1]

uscases=total.cases[usindex]

uscases=as.numeric(as.vector(uscases))

latitude=read.table(file="latitude clusters.txt",header=F,dec=".",sep="\n")

latitude=latitude[,1]

latitude=latitude[usindex]

longitude=read.table(file="longitude clusters.txt",header=F,dec=".")

longitude=longitude[,1]

longitude=longitude[usindex]

#points(longitude,latitude,cex=0.6,col="red")

Y=uscases

####################################

californiaindex=(latitude<=41.21)*(latitude>38.55)*(longitude<=-73.53)*(longitude>-75.35)

Ycalifornia=Y[which(californiaindex==1)]

ntot=length(Ycalifornia)

Zfull=Ycalifornia

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
  
  seuil=which(vecu>=75)[1]
  
  xihat[j]=vecgamma[seuil]
  
  print(j*100/nrep)
  
  j=j+1
}

#save.image(file="robust_nj_13_07.Rdata")