rm(list=ls())

library(maps)

library(evt0)

library(Expectrem)

par(mfrow=c(1,3))

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

californiaindex=(latitude<=46.15)*(latitude>42)*(longitude<=-116.45)*(longitude>-124.30)

Yoregon=Y[which(californiaindex==1)]

points(longitude[which(californiaindex==1)],latitude[which(californiaindex==1)],cex=5*Yoregon/max(Yoregon),col="red")

n=length(Yoregon)

####################################

#################discrete GPD###############

Z=Yoregon

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
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05)
  
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

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="xi",main="SARS-CoV-2 clusters in Oregon",ylim=c(-1,3),xlim=c(0,100),col="red")

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.3,b=0,col="red",lty=5)

abline(v=31,col="blue",lty=5)

tau=0.99

vecquant=1-(1+vecgamma[14]*((1:1000)+1)/vecsigma[14])^(-1/vecgamma[14])

quantau=75+which(vecquant>=tau)[1]

quantau


###############################

n=length(Ztot)

vechill=rep(0,n-1)

vecutot=rep(0,n-1)

k=2

while(k<=n){
  
  vecutot[k-1]=rev(sort(Ztot))[k]
  
  vechill[k-1]=mean(log(Ztot[which(Ztot>=vecutot[k-1])]/vecutot[k-1]))
  
  k=k+1
}

vecgammahill=unique(rev(vechill))

vecuhill=unique(rev(vecutot))


i=1

veck=rep(0,length(vecutot))

while(i<=length(vecutot)){
  
  veck[i]=which(vecutot==vecutot[i])[1]
  
  i=i+1
}

veck=rev(unique(veck))


###################################


rm(list=ls())

library(maps)

library(mapdata)

country=read.table(file="country cluster.txt",sep="\n")

country=country[,1]

sum(country=="Singapore")

usindex=which(country=="Singapore")

total.cases=read.table(file="total cases.txt",header=F,sep="\n")

total.cases=total.cases[,1]

uscases=total.cases[usindex]

uscases=as.numeric(as.vector(uscases))

Z=uscases

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
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05)
  
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

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="xi",main="SARS-CoV-2 clusters in Singapore",ylim=c(-1,3),xlim=c(1,75),col="red")

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.89,b=0,col="red",lty=5)

abline(v=20,col="blue",lty=5)

tau=0.99

vecquant=1-(1+vecgamma[17]*((1:1000)+1)/vecsigma[17])^(-1/vecgamma[17])

quantau=25+which(vecquant>=tau)[1]

quantau


###############################

n=length(Ztot)

vechill=rep(0,n-1)

vecutot=rep(0,n-1)

k=2

while(k<=n){
  
  vecutot[k-1]=rev(sort(Ztot))[k]
  
  vechill[k-1]=mean(log(Ztot[which(Ztot>=vecutot[k-1])]/vecutot[k-1]))
  
  k=k+1
}

vecgammahill=unique(rev(vechill))

vecuhill=unique(rev(vecutot))



rm(list=ls())

library(maps)

library(mapdata)

country=read.table(file="country cluster.txt",sep="\n")

country=country[,1]

sum(country=="South Korea")

usindex=which(country=="South Korea")

total.cases=read.table(file="total cases.txt",header=F,sep="\n")

total.cases=total.cases[,1]

uscases=total.cases[usindex]

uscases=as.numeric(as.vector(uscases))

Z=uscases

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
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05)
  
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

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="xi",main="SARS-CoV-2 clusters in South Korea",ylim=c(-1,3),xlim=c(1,75),col="red")

points(vecgammainfcovid~vecucovid,type="l",lty=3,col="red")

points(vecgammasupcovid~vecucovid,type="l",lty=3,col="red")

abline(a=0.95,b=0,col="red",lty=5)

abline(v=22,col="blue",lty=5)

tau=0.99

vecquant=1-(1+vecgamma[17]*((1:1000)+1)/vecsigma[17])^(-1/vecgamma[17])

quantau=25+which(vecquant>=tau)[1]

quantau


###############################

n=length(Ztot)

vechill=rep(0,n-1)

vecutot=rep(0,n-1)

k=2

while(k<=n){
  
  vecutot[k-1]=rev(sort(Ztot))[k]
  
  vechill[k-1]=mean(log(Ztot[which(Ztot>=vecutot[k-1])]/vecutot[k-1]))
  
  k=k+1
}

vecgammahill=unique(rev(vechill))

vecuhill=unique(rev(vecutot))