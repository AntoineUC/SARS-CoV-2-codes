rm(list=ls())

load(file="traceDatSaved.Rdata")

traceDat <- traceDatSaved
tabTot = table(table(traceDat$id))
tot = c(0,tabTot[1:80],0,0,sum(tabTot[81:length(tabTot)]))
tabTotPos = table(table(traceDat$id[traceDat$cPos==1]))
totPos = rep(0,41); totPos[as.numeric(names(tabTotPos))+1] = tabTotPos
transm = unique(traceDat$id[traceDat$cPos==1])
nontransm = unique(traceDat$id[which((traceDat$id%in%transm)==F)])
p0 = length(nontransm)/(length(transm)+length(nontransm)) # proportion of cases with no secondary infections
tabTotPos = c(length(nontransm),tabTotPos)
names(tabTotPos)[1] = '0'
contactsTot = rep(as.numeric(names(tabTot)),tabTot)
contactsPos = rep(as.numeric(names(tabTotPos)),tabTotPos)
hist(contactsPos)

## Empirical offspring distribution
nsec = contactsPos

Z=nsec

Ztot=Z

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

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
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.975)
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)/sqrt(nk)*qnorm(0.025)
  
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

plot(vecgammacovid~vecucovid,type="l",xlab="Threshold",ylab="gamma",main="SSEs",ylim=c(-0.5,1.5),xlim=c(0,21),col="red")

points(vecgammainfcovid~vecucovid,type="l",col="red",lty=3)

points(vecgammasupcovid~vecucovid,type="l",col="red",lty=3)

abline(a=mean(vecgammacovid[1:11]),b=0,col="red",lty=5)

##########################################

data=data.frame(
  Threshold = vecucovid,
  Tail_Index = vecgammacovid,
  gammasup=vecgammasupcovid,
  gammainf=vecgammainfcovid
)

ggplot(data,aes(x=Threshold, y=Tail_Index)) +
  ylim(c(-1.5, 1)) + xlim(c(0,21)) +
  geom_point(aes(x=Threshold, y=Tail_Index),col="red") +
   geom_line(col="red") +
  geom_point(aes(x=Threshold, y=gammasup),col="red") +
  geom_line(aes(x=Threshold, y=gammasup),col="red",linetype="dashed")+
  geom_point(aes(x=Threshold, y=gammainf),col="red") +
  geom_line(aes(x=Threshold, y=gammainf),col="red",linetype="dashed")+
  geom_abline(intercept = mean(vecgammacovid[1:11]), slope = 0,color="red",linetype="dotted")+
  ylab("Tail  Index") +
  labs(title = "(f) Tail index plot - India data")
