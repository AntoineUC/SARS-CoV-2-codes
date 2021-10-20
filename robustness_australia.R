rm(list=ls())

library(maps)

library(evt0)

library(Expectrem)

dev.off()

#par(mfrow=c(1,1))

country=read.table(file="country cluster.txt",sep="\n")

country=country[,1]

sum(country=="Australia")

usindex=which(country=="Australia")

total.cases=read.table(file="total cases.txt",header=F,sep="\n")

total.cases=total.cases[,1]

uscases=total.cases[usindex]

uscases=as.numeric(as.vector(uscases))

Z=uscases

Z=c(Z,310,301,285,250,240,211,210,197,167,164,160,156,151,147,146,143,141,137,132,132,130,127,127,127,126,118,116,115,113,113,112,111,111,101,96,94,94,93,91,90,89,88,87,86,83,77,76,76,75,75,72,71,71,69,69,68,65,65,65,64,64,64,60,58,55,55,53,53,50,50,49,48,48,48,46,45,45,45,44,44,43,43,41,41,42,40,40,40,39,39,38,37,36,36,35,35,34,34,34,rep(33,6),rep(32,4),rep(31,4),rep(30,2),rep(28,3),rep(27,3),rep(26,2),rep(25,4),rep(24,2),rep(23,4),22,rep(21,6),rep(20,6),rep(19,6),rep(18,5),rep(17,5),rep(16,6),rep(15,4),rep(14,5),rep(13,6),rep(12,8),rep(11,9),rep(10,8),rep(9,10),rep(8,12),rep(7,10),rep(6,14),rep(5,21),rep(4,14),rep(3,14),rep(2,23))

Zfull=Z

ntot=length(Zfull)

xihat=rep(0,10000)

j=1

while(j<=10000){
  
  Ycalifornia=Zfull[-sample(ntot)[1:36]]
  
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

#save.image(file="robust_australia.Rdata")
