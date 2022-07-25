rm(list=ls())

library(ggplot2)

pushViewport(viewport(layout = grid.layout(1,1)))

require(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))

define_region = function(row, col){
  viewport(layout.pos.row=row, layout.pos.col=col)
}

X=read.table(file="sse_korea_2020.txt",header=T)

Y=table(X)

Ztot=Y

Ztot=c(Ztot,rep(0,4558))

x=unique(Ztot)

y=table(factor(Ztot, levels=0:51))

data=as.data.frame(table(factor(Ztot, levels=0:51)))

names(data)=c("Cases","Count")

p=ggplot(data=data, aes(x=Cases, y=log(Count+1))) +
  scale_x_discrete("Cases", breaks=seq(0,51,3))+
  geom_bar(stat="identity", width=0.5)+
  labs(title = "(a) Barplot - Korea data")

print(p,vp=define_region(1,1))

###################################

X=read.table(file="sse_korea_pnas.txt",header=T)

Y=table(X)

Ztot=Y

Ztot=c(Ztot,rep(0,4558))

k=1

#vecu=1:20

vecu=unique(sort(Ztot))

meanexcess=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  meanexcess[k]=mean(Z)
  
  k=k+1
  
}

linreg=lm(meanexcess[1:12]~vecu[1:12])$coefficients

data=data.frame(
  Threshold = vecu,
  Tail_Index = meanexcess
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index)) +
  ylim(c(0, 15)) + xlim(c(0,25)) +
  geom_point(aes(x=Threshold, y=Tail_Index),col="red") +
  #  geom_line(col="red") +
  
  geom_abline(intercept = linreg[1], slope = linreg[2],color="red")+
  ylab("Mean excess value") +
  labs(title = "(b) Mean excess value plots - Korea data") +
  geom_text(x=5, y=12, label=expression(paste(slope,'=0.852',' => ',xi,'=0.460')) ,col="red")

print(p,vp=define_region(1,2))

###########################################L

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

llfunction=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((ddgpd(Z,igamma,sigma)))))
}

dgpd=function(x,igamma,sigma){
  return(1/sigma*(1+igamma*x/sigma)^(-1/igamma-1))
}

llfunctionc=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((dgpd(Zc,igamma,sigma)))))
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

vecuc=seq(min(Ztot),max(Ztot),length=251)

vecgamma=rep(0,length(vecu))

vecgammac=rep(0,length(vecuc))

vecgammah=rep(0,length(vecuc))

vecgammam=rep(0,length(vecuc))

vecgammagh=rep(0,length(vecuc))

vecgammah2=rep(0,length(vecuc))

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
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.975)
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.025)
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.975)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.025)
  
  k=k+1
  
}

k=1

while(k<=length(vecuc)){
  
  Zc=Ztot[which(Ztot>vecuc[k])]-vecuc[k]
  
  Zh=Ztot[which(Ztot>vecuc[k])]/vecuc[k]
  
  thetahatc=optim(fn=llfunctionc,par=c(0.001,0.5))$par
  
  vecgammac[k]=thetahatc[1]
  
  vecgammah[k]=mean(log(Zh))
  
  vecgammah2[k]=mean(log(Zh)^2)
  
  vecgammam[k]=vecgammah[k]+1-0.5*1/(1-vecgammah[k]^2/vecgammah2[k])
  
  k=k+1
  
}

# k=1
# 
# while(k<=length(vecuc)){
#   
#   Zh=Ztot[which(Ztot>vecuc[k])]/vecuc[k]
#   
#   #vecgammagh[k]=sum(log((vecuc[(length(vecuc)-1):(length(vecuc)-k)]*vecgammah[1:k])/(vecuc[(length(vecuc)-k)]*vecgammah[k])))/(sum(Ztot>vecuc[k]))
#   
#   vecgammagh[k]=sum(log((vecuc[(length(vecuc)-1):(length(vecuc)-k)]*vecgammah[1:k])/(vecuc[(length(vecuc)-k)]*vecgammah[k])))/(sum(Ztot>vecuc[k]))
#   
#   k=k+1
#   
# }

vecgamma

vecgammacovid=vecgamma

vecgammainfcovid=vecgammainf

vecgammasupcovid=vecgammasup

vecsigmacovid=vecsigma

vecsigmasupcovid=vecsigmasup

vecsigmainfcovid=vecsigmainf

vecucovid=vecu

data=data.frame(
  Threshold = vecucovid,
  Tail_Index = vecgammacovid
)

data2=data.frame(
  Thresholdc=vecuc,
  GPDc=vecgammac,
  Hillc=vecgammah,
  Momentc=vecgammam
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index)) +
  ylim(c(-1.5, 2.25)) + xlim(c(0,21)) +
  geom_point(aes(x=Threshold, y=Tail_Index),col="red") +
  geom_line(col="red") +
  #geom_point(aes(x=Threshold, y=gammasup),col="red") +
  #geom_line(aes(x=Threshold, y=gammasup),col="red",linetype="dashed")+
  #geom_point(aes(x=Threshold, y=gammainf),col="red") +
  #geom_line(aes(x=Threshold, y=gammainf),col="red",linetype="dashed")+
  #  geom_abline(intercept = 0.209, slope = 0,color="red",linetype="dashed",lwd=1.5)+
  geom_line(col="red",lwd=2) +
  geom_point(data=data2,aes(x=Thresholdc, y=Hillc),col="black") +geom_line(data=data2,aes(x=Thresholdc, y=Hillc),col="black") +
  geom_point(data=data2,aes(x=Thresholdc, y=GPDc),col="blue") +geom_line(data=data2,aes(x=Thresholdc, y=GPDc),col="blue") +
  #  geom_point(data=data2,aes(x=Thresholdc, y=Momentc),col="forestgreen") +geom_line(data=data2,aes(x=Thresholdc, y=Momentc),col="forestgreen") +
  ylab("Tail  Index") +
  labs(title = "(c) Three tail index estimates - Korea data")

print(p,vp=define_region(1,1))
# 
# plot(vecgamma~vecu,type="l",col="red",xlim=c(0,21),ylim=c(-1,1))
# 
# points(vecgammac~vecuc,type="l",col="blue")
# 
# points(vecgammah~vecuc,type="l",col="black")

#####################Bottom panels########################

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

llfunction=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((ddgpd(Z,igamma,sigma)))))
}

dgpd=function(x,igamma,sigma){
  return(1/sigma*(1+igamma*x/sigma)^(-1/igamma-1))
}

llfunctionc=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((dgpd(Zc,igamma,sigma)))))
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

vecuc=seq(min(Ztot),max(Ztot),length=196)

vecgamma=rep(0,length(vecu))

vecgammac=rep(0,length(vecuc))

vecgammah=rep(0,length(vecuc))

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
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.975)
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.025)
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.975)
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)/sqrt(nk)*qnorm(0.025)
  
  k=k+1
  
}

k=1

while(k<=length(vecuc)){
  
  Zc=Ztot[which(Ztot>vecuc[k])]-vecuc[k]
  
  Zh=Ztot[which(Ztot>vecuc[k])]/vecuc[k]
  
  thetahatc=optim(fn=llfunctionc,par=c(0.001,0.5))$par
  
  vecgammac[k]=thetahatc[1]
  
  vecgammah[k]=mean(log(Zh))
  
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

data=data.frame(
  Threshold = vecucovid,
  Tail_Index = vecgammacovid,
  gammasup=vecgammasupcovid,
  gammainf=vecgammainfcovid
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index)) +
  ylim(c(-1.5, 1)) + xlim(c(0,21)) +
  geom_point(aes(x=Threshold, y=Tail_Index),col="red") +
  geom_line(col="red") +
  geom_point(aes(x=Threshold, y=gammasup),col="red") +
  geom_line(aes(x=Threshold, y=gammasup),col="red",linetype="dashed")+
  geom_point(aes(x=Threshold, y=gammainf),col="red") +
  geom_line(aes(x=Threshold, y=gammainf),col="red",linetype="dashed")+
  #  geom_abline(intercept = 0.209, slope = 0,color="red",linetype="dotted")+
  #  ylab("Tail  Index") +
  labs(title = "(d) Tail index plot D-GPD - Korea data")

print(p,vp=define_region(1,2))