rm(list=ls())

library(maps)

library(evt0)

library(Expectrem)

dev.off()

#par(mfrow=c(1,1))

country=read.table(file="country cluster.txt",sep="\n")

country=country[,1]

sum(country=="Canada")

usindex=which(country=="Canada")

total.cases=read.table(file="total cases.txt",header=F,sep="\n")

total.cases=total.cases[,1]

uscases=total.cases[usindex]

uscases=as.numeric(as.vector(uscases))

Z=uscases

Zfull=Z

ntot=length(Zfull)

xihat=rep(0,10000)

j=1

while(j<=10000){
  
  Ycalifornia=Zfull[-sample(ntot)[1:10]]
  
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
  
  seuil=which(vecu>=25)[1]
  
  xihat[j]=vecgamma[seuil]
  
  print(j*100/10000)
  
  j=j+1
}

#save.image(file="robust_canada.Rdata")
