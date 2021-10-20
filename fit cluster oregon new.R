rm(list=ls())

library(maps)

library(ggplot2)

library(mapdata)

par(mfrow=c(1,1))

setwd("C:/Users/antoi/Dropbox/Post doc ENSAI/PNAS")

Z=c(7,5,9,46,19,7,6,9,4,6,6,3,10,8,1,19,5,4,17,13,6,4,3,43,17,15,11,2,4,9,6,47,6,6,9,6,6,3,4,134,3,8,7,2,4,3,4,29
    ,3,78,3,51,3,13,31,4,8,37,5,8,5,6,2,5,2,4,4,25,4,51,7,33,141,5,6,9,14,8,3,16,44,4,4,8,5,3,3,14,3,23,24,8,34,17
    ,3,33,1,3,29,33,10,5,6,62,4,5,3,3,37,4,26,3,14,19,15,12,4,8,68,37,3,12,5,31,13,13,22,37,4,4,7,6,15,67,4,3,3
    ,9,19,4,3,5,35,8,76,15,33,18,22,8,56,4,34,3,3,8,3,40,9,2,6,37,3,7,6,15,30,14,48,17,4,18,21,121,3,8,8,3,12
    ,12,41,4,4,3,4,6,2,4,10,8,5,7,3,3,50,31,53,1,11,14,65,4,6,10,93,23,3,13,5,5,5,27,7,16,6,7,41,4,47,5,25,7,13,1
    ,40,14,30,23,6,41,34,6,8,12,3,5,15,2,31,6,5,23,4,51,60,4,79,67,19,3,4,5,1,5,92,15,4,3,38,17,3,9,13,13,13,4,21
    ,11,3,8,6,4,31,4,11,112,3,9,12,12,20,5,3,18,12,58,4,5,6,26,24,116,14,50,4,103,6,11,23,2,3,21,15,7,15,4,40,9,55
    ,28,3,9,9,4,5,3,71,2,26,3,3,30,2,51,9,9,7,8,11,47,27,4,6,13,7,3,7,58,3,6,4,6,6,14,10,38,3,83,20,3,53,19,10,58
    ,4,3,5,24,7,29,16,16,22,4,46,10,3,5,111,8,6,5,113,49,3,24,7,6,6,5,23,6,5,3,31,3,4,71,37,12,3,11,4,104,1,5,112
    ,6,6,3,7,29,3,2,60,27,5,11,5,3,11,7,11,8,3,15,5,4,5,7,22,11,21,3,4,34,22,25,30,66,24,9,34,3,1,5,3,9,72,13,14,27
    ,3,8,7,9,3,31,10,39,3,13,3,52,21,11,4,3,56,4,15,3,3,29,5,26,18,41,48,7,4,24,2,8,4,18,2,4,7,9,34,60,3,39,9,27,4
    ,3,40,30,10,30,3,4,93,3,57,3,88,10,5,13,55,3,11,6,9,4,3,17,6,7,39,7,34,9,3,3,3,4,4,3,3,6,5,13,4,4,49,2,7,32,6
    ,6,4,34,11,3,4,7,31,97,6,21,21,2,58,14,5,29,3,5,7,3,4,10,8,26,6,13,7,3,3,5,62,2,10,4,38,5,70,4,57,5,19,8,4,14
    ,4,7,3,5,50,5,28,4,13,3,48,13,15,7,6,3,3,22,47,3,32,24,26,8,32,3,92,45,23,27,53,37,11,21,23,30,6,11,17,4,6,3
    ,5,5,5,18,32,15,26,8,28,4,3,11,6,3,7,9,23,9,10,4,6,14,4,5,5,22,41,19,11,37,3,2,3,54,6,25,10,53,23,5,4,9,41,20
    ,7,27,18,113,24,32,7,6,6,4,12,11,32,20,8,4,4,9,4,3,14,8,22,18,10,3,10,3,7,21,14,3,3,8,6,4,16,41,3,4,8,8,4,17
    ,3,3,6,3,3,6,6,3,3,22,53,4,73,10,8,26,27,3,3,3,5,11,7,6,11,4) #resolved care facilities outbreaks

Z=c(Z,32,22,16,10,10,10,9,9,8,7,7,7,rep(6,5),rep(5,10)) #resolved workplace outbreaks 01/09

Z=c(Z,34,14,13,11,9,8,6,6,6,6,5,5) #resolved workplace outbreaks 25/08

Z=c(Z,22,11,8,7,6,6,6,5,5) #resolved workplace outbreaks 18/08

Z=c(Z,639,337,31,11,8,6) #resolved workplace outbreaks 11/08

Z=c(Z,38,27,11,9,7,7,6) #resolved workplace outbreaks 04/08

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

n=length(Z)

Ztot=Z

########################################

vecu=unique(sort(Z))

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}

plot(y=meanexcess,x=vecu,pch=1,col="red",ylim=c(0,300),xlim=c(0,250),ylab="Mean excess value",xlab="Threshold (u)",main="SARS-CoV-2 clusters in Hong-Kong")

linreg=lm(meanexcess[1:25]~vecu[1:25])$coefficients

abline(a=linreg[1],b=linreg[2],col="red")

linreg[2]/(1+linreg[2])

text(x=450,y=250,labels="a=8.949 => gamma=0.899",col="red")

###################################

llfunction=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((ddgpd(Z,igamma,sigma)))))
}

llfunction2=function(igamma,sigma){
  return(-sum(log(((1+igamma*Z/sigma)^(-1/igamma)-(1+igamma*(Z+1)/sigma)^(-1/igamma)))))
}

dllfunction2gamma=function(igamma,sigma){
  return((llfunction2(igamma+0.001,sigma)-llfunction2(igamma,sigma))/0.001)
}

dllfunction2gammasigma=function(igamma,sigma){
  return((dllfunction2gamma(igamma,sigma+0.001)-dllfunction2gamma(igamma,sigma))/0.001)
}

dllfunction2gamma2=function(igamma,sigma){
  return((dllfunction2gamma(igamma+0.001,sigma)-dllfunction2gamma(igamma,sigma))/0.001)
}

dllfunction2sigma=function(igamma,sigma){
  return((llfunction2(igamma,sigma+0.001)-llfunction2(igamma,sigma))/0.001)
}

dllfunction2sigma2=function(igamma,sigma){
  return((dllfunction2sigma(igamma,sigma+0.001)-dllfunction2sigma(igamma,sigma))/0.001)
}

k=1

#vecu=1:20

vecu=unique(sort(Ztot))

vecgamma=rep(0,length(vecu))

vecgammainf=rep(0,length(vecu))

vecgammasup=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecsigmainf=rep(0,length(vecu))

vecsigmasup=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  nk=length(Z)
  
  thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par
  
  vecgamma[k]=thetahat[1]
  
  vecsigma[k]=thetahat[2]
  
  I11=dllfunction2gamma2(thetahat[1],thetahat[2])
  
  I12=dllfunction2gammasigma(thetahat[1],thetahat[2])
  
  I21=I12
  
  I22=dllfunction2sigma2(thetahat[1],thetahat[2])
  
  matI=matrix(c(I11,I12,I21,I22),2,2)
  
  varhat=solve(matI)[1,1]
  
  varhatsigma=solve(matI)[2,2]
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.95)
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.05)
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.975)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.025)
  
  k=k+1
  
}

vecgamma

vecgammacovid=vecgamma

vecgammainfcovid=vecgammainf

vecgammasupcovid=vecgammasup

vecsigmacovid=vecsigma

vecsigmasupcovid=vecsigmasup

vecsigmainfcovid=vecsigmainf

vecucovid=vecu

dev.off()

par(mfrow=c(1,1))

hist(Ztot,breaks=50,freq=T,main="",ylim=c(0,150),xlim=c(1,650),xlab="",ylab="")

title("SARS-CoV-2 clusters in Colorado",line=2.5)

par(new=T)

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="Count",col="red",ylim=c(-1,2),xlim=c(0,50),axes=F)

axis(side=4, at = (-2:4)/2)

mtext("xi", side=4, line=-1.5)

axis(side=3, at = (0:10)*5)

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.21,b=0,col="red",lty=5)

abline(v=15,col="red",lty=5)

tau=0.99

vecquant=((1+vecgamma[15]*((0:10000))/vecsigma[15])^(-1/vecgamma[15]))*mean(Ztot>=vecu[15])

quantau=15-1+which(vecquant<=1-tau)[1]

quantau