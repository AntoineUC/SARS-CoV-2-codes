rm(list=ls())

dev.off()

par(mfrow=c(1,1))

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

####################################

#n=1000

#Z=sample(x = 0:100, n, replace = T, prob = ddgpd(0:100,0.5,1))

####################################

Z=c(580,68,39,275,46,61,40,38,73,64,134,21,1,1,38,50,9,26,30,197,57,12,22,113,26,2,7,19,68,39,46,33,87,85,39,24,19,49,55,31,48,132,63,27,17,54,46,12,13,22,9,44,36,49,47,61,17,12,23,12,4,13,13,14,13,4,16,6,6,6,17,11,66,15,4,3,3,2,2,2,13,18,33,17,31,20,18,15,38,7,13,21,3,18,22,24,84,39,22,20,18,67,66,27,18,34,41,17,20,29,35,17,12)

Z=Z[,1]

#Z=c(Z,57)

Ztot=Z

########################################

vecu=unique(sort(Z))

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}

plot(y=meanexcess,x=vecu,pch=1,col="red",ylim=c(0,200),xlim=c(0,300),ylab="Mean excess value",xlab="Threshold (u)",main="SARS-CoV-2 clusters in Kerala")

linreg=lm(meanexcess[5:52]~vecu[5:52])$coefficients

abline(a=linreg[1],b=linreg[2],col="red")

linreg[2]/(1+linreg[2])

text(x=170,y=30,labels="a=0.672 => gamma=0.402",col="red")


###################################

# Z=read.table(file="korea_pnas.txt",header=F)
# 
# Z=Z[,1]
# 
# n=length(Z)
# 
# Ztot=Z

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

vecmoment=rep(0,length(vecu))

vecsigmainf=rep(0,length(vecu))

vecsigmasup=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  nk=length(Z)
  
  thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par
  
  vecmoment[k]=0.5*(1-(mean(Z)/sd(Z))^2)
  
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

hist(Ztot,breaks=25,freq=T,main="",ylim=c(0,40),xlim=c(1,600),xlab="",ylab="")

title("SARS-CoV-2 clusters in Kerala",line=2.5)

par(new=T)

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="Count",col="red",ylim=c(-1,2),xlim=c(0,120),axes=F)

axis(side=4, at = (-2:4)/2)

mtext("xi", side=4, line=-1.5)

axis(side=3, at = (0:10)*12)

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

#points(vecmoment~vecucovid,type="l",col="green")

abline(a=0.45,b=0,col="red",lty=5)

abline(v=27,col="red",lty=5)

tau=0.95

vecquant=((1+vecgamma[19]*((0:10000))/vecsigma[19])^(-1/vecgamma[19]))*mean(Ztot>=vecu[19])

quantau=22-1+which(vecquant<=1-tau)[1]

quantau
