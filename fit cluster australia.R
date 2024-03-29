rm(list=ls())

library(maps)

library(ggplot2)

library(mapdata)

par(mfrow=c(1,1))

load(file="clusters.Rdata")

country=clusters[,2]

sum(country=="Australia")

usindex=which(country=="Australia")

uscases=clusters[usindex,1]

uscases=as.numeric(as.vector(uscases))

Z=uscases

Z=c(Z,310,301,285,250,240,211,210,197,167,164,160,156,151,147,146,143,141,137,132,132,130,127,127,127,126,118,116,115,113,113,112,111,111,101,96,94,94,93,91,90,89,88,87,86,83,77,76,76,75,75,72,71,71,69,69,68,65,65,65,64,64,64,60,58,55,55,53,53,50,50,49,48,48,48,46,45,45,45,44,44,43,43,41,41,42,40,40,40,39,39,38,37,36,36,35,35,34,34,34,rep(33,6),rep(32,4),rep(31,4),rep(30,2),rep(28,3),rep(27,3),rep(26,2),rep(25,4),rep(24,2),rep(23,4),22,rep(21,6),rep(20,6),rep(19,6),rep(18,5),rep(17,5),rep(16,6),rep(15,4),rep(14,5),rep(13,6),rep(12,8),rep(11,9),rep(10,8),rep(9,10),rep(8,12),rep(7,10),rep(6,14),rep(5,21),rep(4,14),rep(3,14),rep(2,23))

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

plot(y=meanexcess,x=vecu,pch=1,col="red",ylim=c(0,300),xlim=c(0,250),ylab="Mean excess value",xlab="Threshold (u)",main="SARS-CoV-2 clusters in Australia")

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

hist(Ztot,breaks=25,freq=T,main="",ylim=c(0,200),xlim=c(1,700),xlab="",ylab="")

title("SARS-CoV-2 clusters in Australia",line=2.5)

par(new=T)

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="Count",col="red",ylim=c(-1,2),xlim=c(0,50),axes=F)

axis(side=4, at = (-2:4)/2)

mtext("xi", side=4, line=-1.5)

axis(side=3, at = (0:10)*5)

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.28,b=0,col="red",lty=5)

abline(v=25,col="red",lty=5)

tau=0.95

vecquant=((1+vecgamma[24]*((0:10000))/vecsigma[24])^(-1/vecgamma[24]))*mean(Ztot>=vecu[24])

quantau=25-1+which(vecquant<=1-tau)[1]

quantau
