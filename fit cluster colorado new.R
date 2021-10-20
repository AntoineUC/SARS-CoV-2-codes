rm(list=ls())

library(maps)

library(ggplot2)

library(mapdata)

par(mfrow=c(1,1))

setwd("C:/Users/antoi/Dropbox/Post doc ENSAI/PNAS")

Z=read.table(file="colorado_new.txt")

Z=Z[,1]

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

hist(Ztot,breaks=50,freq=T,main="",ylim=c(0,150),xlim=c(1,1750),xlab="",ylab="")

title("SARS-CoV-2 clusters in Colorado",line=2.5)

par(new=T)

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="Count",col="red",ylim=c(-1,2),xlim=c(0,50),axes=F)

axis(side=4, at = (-2:4)/2)

mtext("xi", side=4, line=-1.5)

axis(side=3, at = (0:10)*5)

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.53,b=0,col="red",lty=5)

abline(v=28,col="red",lty=5)

tau=0.99

vecquant=((1+vecgamma[26]*((0:10000))/vecsigma[26])^(-1/vecgamma[26]))*mean(Ztot>=vecu[26])

quantau=27-1+which(vecquant<=1-tau)[1]

quantau