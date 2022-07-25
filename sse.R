rm(list=ls())

library(ggplot2)

require(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,2)))

define_region = function(row, col){
  viewport(layout.pos.row=row, layout.pos.col=col)
}

############################################

#data

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

#####################################

data <- data.frame(
  Disease = c( rep("SARS-CoV ", 15), rep("SARS-CoV-2", 45) ),
  Cases = c( Z[46:60], Z[1:45] )
)

p=ggplot(data, aes(x=Cases,fill=Disease)) +
  geom_histogram(alpha=0.5,position="identity")


print(p+labs(title = "(a) Histograms"),vp=define_region(1,1))

###########################################


Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Zcov=Z[1:45] #SARS-CoV-2

Zsars=Z[46:60] #SARS-CoV

f=function(theta){
  return(-sum(log((dnbinom(Zcov,size=theta[1],mu=theta[2]))/(1-pnbinom(6,size=theta[1],mu=theta[2])))))
}

thetahatcov=optim(fn=f,par=c(100,100))$par

cumsum(dnbinom(7:100,size=thetahatcov[1],mu=thetahatcov[2])/(1-pnbinom(6,size=thetahatcov[1],mu=thetahatcov[2])))

#######################################


f=function(theta){
  return(-sum(log((dnbinom(Zsars,size=theta[1],mu=theta[2]))/(1-pnbinom(6,size=theta[1],mu=theta[2])))))
}

thetahatsars=optim(fn=f,par=c(100,100))$par

cumsum(dnbinom(7:100,size=thetahatsars[1],mu=thetahatsars[2])/(1-pnbinom(6,size=thetahatsars[1],mu=thetahatsars[2])))

data <- data.frame(
  Disease = c( rep("SARS-CoV-2", 34), rep("SARS-CoV ", 34) ),
  probs = c((dnbinom(7:40,size=thetahatcov[1],mu=thetahatcov[2])/(1-pnbinom(6,size=thetahatcov[1],mu=thetahatcov[2]))),(dnbinom(7:40,size=thetahatsars[1],mu=thetahatsars[2])/(1-pnbinom(6,size=thetahatsars[1],mu=thetahatsars[2])))),
  xaxis = c((7:40),(7:40))
)

p=ggplot(data,aes(x = xaxis, y = probs, fill = Disease)) +
  geom_col(alpha=0.6,position="identity") +
  labs(title = "(b)  Negative  binomial  fits",
       x = "",
       y = "Probability") 

print(p,vp=define_region(1,2))


###########################################



Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Ztot=Z

Z=Z[46:60] 

########################################

Z=Ztot[1:45] #SARS-CoV-2

vecucov2=sort(Z)

meanexcesscov2=rep(0,length(vecucov2))

i=1

while(i<=length(vecucov2)){
  
  meanexcesscov2[i]=mean(Z[which(Z>=vecucov2[i])]-vecucov2[i])
  
  i=i+1
}

Z=Ztot[46:60] #SARS-CoV

vecu=sort(Z)

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}

#####################################

data <- data.frame(
  Disease = c( rep("SARS-CoV ", 15), rep("SARS-CoV-2", 45) ),
  Excess = c( meanexcess, meanexcesscov2 ),
  Threshold = c( vecu , vecucov2 )
)

linreg=lm(meanexcess[1:13]~vecu[1:13])$coefficients

linreg2=lm(meanexcesscov2[1:43]~vecucov2[1:43])$coefficients

p=ggplot(data, aes(x=Threshold, y=Excess,col=Disease)) +
  scale_color_manual(values = c("red", "blue")) +
  geom_point() + ylim(c(0, 100)) + xlim(c(0,50))+ 
  geom_abline(intercept = linreg[1], slope = linreg[2],color="blue") +
  geom_abline(intercept = linreg2[1], slope = linreg2[2],color="red") +
  annotate(geom="text",x=40, y=50, label=expression(paste(slope,'=1.43',' => ',xi,'=0.589')),color="blue") +
  annotate(geom="text",x=40, y=5, label=expression(paste(slope,'=0.256',' => ',xi,'=0.204')),color="red") +
  ylab("Mean  Excess  Value")

print(p+labs(title = "(d) Mean Excess Value plots")
      ,vp=define_region(2,2))

##############################################



Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z60=Z

Ztot=Z

Z=Z[46:60]

########################################

Z=Ztot[1:45] #SARS-CoV-2

vecu=sort(Z)

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}


linreg=lm(meanexcess[1:43]~vecu[1:43])$coefficients


linreg[2]/(1+linreg[2]) #approximation of xi

Z=Ztot[46:60] #SARS-CoV

vecu=sort(Z)

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}

linreg=lm(meanexcess[1:13]~vecu[1:13])$coefficients

abline(a=linreg[1],b=linreg[2],col="blue")

linreg[2]/(1+linreg[2]) #approximation of the xi

#####################################


ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
} #mass function


Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[1:45] #SARS-CoV-2

Ztot=Z

################derivatives approximations###################

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

vecgammahillcovid=rep(0,length(vecu))

vecgammasuphillcovid=rep(0,length(vecu))

vecgammainfhillcovid=rep(0,length(vecu))

vecgammainf=rep(0,length(vecu))

vecgammasup=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecsigmainf=rep(0,length(vecu))

vecsigmasup=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  nk=length(Z)
  
  thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par #D-GPD fit
  
  vecgammahillcovid[k]=mean(log(Ztot[which(Ztot>=rev(vecu)[k])]/rev(vecu)[k])) #hill estimator
  
  vecgammasuphillcovid[k]=vecgammahillcovid[k]*(1+qnorm(0.95)/sqrt(nk))
  
  vecgammainfhillcovid[k]=vecgammahillcovid[k]*(1+qnorm(0.05)/sqrt(nk))
  
  vecgamma[k]=thetahat[1]
  
  vecsigma[k]=thetahat[2]
  
  I11=dllfunction2gamma2(thetahat[1],thetahat[2])
  
  I12=dllfunction2gammasigma(thetahat[1],thetahat[2])
  
  I21=I12
  
  I22=dllfunction2sigma2(thetahat[1],thetahat[2])
  
  matI=matrix(c(I11,I12,I21,I22),2,2)
  
  varhat=solve(matI)[1,1]
  
  varhatsigma=solve(matI)[2,2]
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.95) #upper bound for xi
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.05) #lower bound for xi
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95) #upper bound for sigma
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05) #lower bound for sigma
  
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

####################################

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[46:60] #SARS-CoV

Ztot=Z

###################################

llfunction=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((ddgpd(Z,igamma,sigma)))))
} #log-likelihood function

k=1

vecu=unique(sort(Ztot))

vecu=c(6,vecu)

vecgamma=rep(0,length(vecu))

vecgammainf=rep(0,length(vecu))

vecgammasup=rep(0,length(vecu))

vecgammahillsars=rep(0,length(vecu))

vecgammainfhillsars=rep(0,length(vecu))

vecgammasuphillsars=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecsigmasup=rep(0,length(vecu))

vecsigmainf=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  nk=length(Z)
  
  thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par
  
  vecgammahillsars[k]=mean(log(Ztot[which(Ztot>=rev(vecu)[k])]/rev(vecu)[k])) #hill estimator
  
  vecgammasuphillsars[k]=vecgammahillsars[k]*(1+qnorm(0.95)/sqrt(nk))
  
  vecgammainfhillsars[k]=vecgammahillsars[k]*(1+qnorm(0.05)/sqrt(nk))
  
  vecgamma[k]=thetahat[1]
  
  vecsigma[k]=thetahat[2]
  
  I11=dllfunction2gamma2(thetahat[1],thetahat[2])
  
  I12=dllfunction2gammasigma(thetahat[1],thetahat[2])
  
  I21=I12
  
  I22=dllfunction2sigma2(thetahat[1],thetahat[2])
  
  matI=matrix(c(I11,I12,I21,I22),2,2)
  
  varhat=solve(matI)[1,1]
  
  varhatsigma=solve(matI)[2,2]
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.95) #upper bound for xi
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.05) #lower bound for xi
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95) #upper bound for sigma
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05) #lower bound for sigma
  
  k=k+1
  
}

vecgamma

vecgammasars=vecgamma

vecgammainfsars=vecgammainf

vecgammasupsars=vecgammasup

vecsigmasars=vecsigma

vecsigmainfsars=vecsigmainf

vecsigmasupsars=vecsigmasup

vecusars=vecu

########################################

Z60=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

n=length(Z60)

vechill=rep(0,n-1)

vecu60=rep(0,n-1)

k=2

while(k<=n){
  
  vecu60[k-1]=rev(sort(Z60))[k]
  
  vechill[k-1]=mean(log(Z60[which(Z60>=vecu60[k-1])]/vecu60[k-1])) #hill estimator
  
  k=k+1
}

vecgammahill=unique(rev(vechill))

vecuhill=unique(rev(vecu60))

data <- data.frame(
  Disease = c( rep("SARS-CoV ", 13), rep("SARS-CoV-2", 17), rep("Pooled data",59) ),
  Tail_Index = c( vecgammahillsars , vecgammahillcovid , vechill ),
  Threshold = c( vecusars , vecucovid , vecu60 ),
  CIdown=c( vecgammainfhillsars , vecgammainfhillcovid , rep(-10,59) ),
  CIup=c( vecgammasuphillsars , vecgammasuphillcovid , rep(-10,59) )
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index,col=Disease)) +
  scale_color_manual(values = c("black", "red", "blue")) +
  geom_point(aes(x=Threshold, y=Tail_Index,col=Disease)) + geom_line() +
  geom_point(aes(x=Threshold, y=CIdown,col=Disease))  + geom_line(aes(x=Threshold, y=CIdown,col=Disease),linetype="dashed") +
  geom_point(aes(x=Threshold, y=CIup,col=Disease))  + geom_line(aes(x=Threshold, y=CIup,col=Disease),linetype="dashed") +
  ylim(c(-0.5, 1.5)) + xlim(c(6,25)) +
  geom_line() + 
  #geom_abline(intercept = 0.25, slope = 0,color="red",linetype="dotted") +
  #geom_abline(intercept = 0.8, slope = 0,color="blue",linetype="dotted") +
  #geom_abline(intercept = 0.6, slope = 0,color="black",linetype="dotted") + 
  #annotate(geom="text",x=22, y=0.8, label=expression(paste(xi,'=0.8')),color="blue") +
  #annotate(geom="text",x=22, y=0.6, label=expression(paste(xi,'=0.6')),color="black") +
  #annotate(geom="text",x=22, y=0.25, label=expression(paste(xi,'=0.25')),color="red") +
  ylab("Tail  Index") 


print(p+labs(title = "(c) Tail index plots")#+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      #                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      ,vp=define_region(2,1))

###########################################

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z60=Z

Ztot=Z

Z=Z[46:60]

########################################

Z=Ztot[1:45] #SARS-CoV-2

vecu=sort(Z)

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}


linreg=lm(meanexcess[1:43]~vecu[1:43])$coefficients


linreg[2]/(1+linreg[2]) #approximation of xi

Z=Ztot[46:60] #SARS-CoV

vecu=sort(Z)

meanexcess=rep(0,length(vecu))

i=1

while(i<=length(vecu)){
  
  meanexcess[i]=mean(Z[which(Z>=vecu[i])]-vecu[i])
  
  i=i+1
}

linreg=lm(meanexcess[1:13]~vecu[1:13])$coefficients

abline(a=linreg[1],b=linreg[2],col="blue")

linreg[2]/(1+linreg[2]) #approximation of the xi

#####################################

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
} #mass function


Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[1:45] #SARS-CoV-2

Ztot=Z

################derivatives approximations###################

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
  
  thetahat=optim(fn=llfunction,par=c(0.25,0.5))$par #D-GPD fit
  
  vecgamma[k]=thetahat[1]
  
  vecsigma[k]=thetahat[2]
  
  I11=dllfunction2gamma2(thetahat[1],thetahat[2])
  
  I12=dllfunction2gammasigma(thetahat[1],thetahat[2])
  
  I21=I12
  
  I22=dllfunction2sigma2(thetahat[1],thetahat[2])
  
  matI=matrix(c(I11,I12,I21,I22),2,2)
  
  varhat=solve(matI)[1,1]
  
  varhatsigma=solve(matI)[2,2]
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.95) #upper bound for xi
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.05) #lower bound for xi
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95) #upper bound for sigma
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05) #lower bound for sigma
  
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

####################################

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[46:60] #SARS-CoV

Ztot=Z

###################################

llfunction=function(theta){
  igamma=theta[1]
  sigma=theta[2]
  return(-sum(log((ddgpd(Z,igamma,sigma)))))
} #log-likelihood function

k=1

vecu=unique(sort(Ztot))

vecu=c(6,vecu)

vecgamma=rep(0,length(vecu))

vecgammainf=rep(0,length(vecu))

vecgammasup=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecsigmasup=rep(0,length(vecu))

vecsigmainf=rep(0,length(vecu))

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
  
  vecgammasup[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.95) #upper bound for xi
  
  vecgammainf[k]=vecgamma[k]+sqrt(varhat)*qnorm(0.05) #lower bound for xi
  
  vecsigmasup[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.95) #upper bound for sigma
  
  vecsigmainf[k]=vecsigma[k]+sqrt(varhatsigma)*qnorm(0.05) #lower bound for sigma
  
  k=k+1
  
}

vecgamma

vecgammasars=vecgamma

vecgammainfsars=vecgammainf

vecgammasupsars=vecgammasup

vecsigmasars=vecsigma

vecsigmainfsars=vecsigmainf

vecsigmasupsars=vecsigmasup

vecusars=vecu

########################################

n=length(Z60)

vechill=rep(0,n-1)

vecu60=rep(0,n-1)

k=2

while(k<=n){
  
  vecu60[k-1]=rev(sort(Z60))[k]
  
  vechill[k-1]=mean(log(Z60[which(Z60>=vecu60[k-1])]/vecu60[k-1])) #hill estimator
  
  k=k+1
}

vecgammahill=unique(rev(vechill))

vecuhill=unique(rev(vecu60))

data <- data.frame(
  Disease = c( rep("SARS-CoV ", 13), rep("SARS-CoV-2", 17), rep("Pooled data",22) ),
  Tail_Index = c( vecgammasars , vecgammacovid , vecgammahill ),
  Threshold = c( vecusars , vecucovid , vecuhill ),
  CIdown=c( vecgammainfsars , vecgammainfcovid , rep(-10,22) ),
  CIup=c( vecgammasupsars , vecgammasupcovid , rep(-10,22) )
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index,col=Disease)) +
  scale_color_manual(values = c("black", "red", "blue")) +
  geom_point(aes(x=Threshold, y=Tail_Index,col=Disease)) + geom_line() +
  geom_point(aes(x=Threshold, y=CIdown,col=Disease))  + geom_line(aes(x=Threshold, y=CIdown,col=Disease),linetype="dashed") +
  geom_point(aes(x=Threshold, y=CIup,col=Disease))  + geom_line(aes(x=Threshold, y=CIup,col=Disease),linetype="dashed") +
  ylim(c(-1, 3)) + xlim(c(5,22)) +
  geom_line() + 
  #geom_abline(intercept = 0.25, slope = 0,color="red",linetype="dotted") +
  #geom_abline(intercept = 0.8, slope = 0,color="blue",linetype="dotted") +
  #geom_abline(intercept = 0.6, slope = 0,color="black",linetype="dotted") + 
  #annotate(geom="text",x=22, y=0.8, label=expression(paste(xi,'=0.8')),color="blue") +
  #annotate(geom="text",x=22, y=0.6, label=expression(paste(xi,'=0.6')),color="black") +
  #annotate(geom="text",x=22, y=0.25, label=expression(paste(xi,'=0.25')),color="red") +
  ylab("Tail  Index") 


print(p+labs(title = "(e) Tail index plots D-GPD")#+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      #                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      ,vp=define_region(3,1))

#################################################################


sigma=1

igamma=0.5

ddgpd=function(k,igamma,sigma){
  return((1+igamma*k/sigma)^(-1/igamma)-(1+igamma*(k+1)/sigma)^(-1/igamma))
}

####################################

#n=1000

#Z=sample(x = 0:100, n, replace = T, prob = ddgpd(0:100,0.5,1))

####################################

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[1:45]

Ztot=Z

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

k=1

#vecu=1:20

vecu=unique(sort(Ztot))

vecgamma=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecmass=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  vecgamma[k]=optim(fn=llfunction,par=c(0.25,0.5))$par[1]
  
  vecsigma[k]=optim(fn=llfunction,par=c(0.25,0.5))$par[2]
  
  vecmass[k]=ddgpd(vecu[k],vecgamma[k],vecsigma[k])
  
  k=k+1
  
}

vecgamma

vecgammacovid=vecgamma

vecsigmacovid=vecsigma

vecmasscovid=vecmass

vecucovid=vecu

####################################

Z=c(10,9,6.5,11,11,8,7,7,25,8,7,9,12,27,21,17,15,11,10,13,7,17.5,52,7,8,7,18,24,11,51,21,15,7,9,17,12,24,8,11,10,10,8,8,11,15,21,19,13,13,33,112,21,23,23,40,187,15,10,8,12)

Z=round(Z)

Z=Z[46:60]

Ztot=Z

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

k=1

#vecu=1:20

vecu=unique(sort(Ztot))

vecu=c(6,vecu)

vecgamma=rep(0,length(vecu))

vecsigma=rep(0,length(vecu))

vecmass=rep(0,length(vecu))

while(k<=length(vecu)){
  
  Z=Ztot[which(Ztot>=vecu[k])]-vecu[k]
  
  vecgamma[k]=optim(fn=llfunction,par=c(0.25,0.5))$par[1]
  
  vecsigma[k]=optim(fn=llfunction,par=c(0.25,0.5))$par[2]
  
  vecmass[k]=ddgpd(vecu[k],vecgamma[k],vecsigma[k])
  
  k=k+1
  
}

vecgamma

vecgammasars=vecgamma

vecsigmasars=vecsigma

vecmasssars=vecmass

vecusars=vecu

ycovid=log(ddgpd(0:100,vecgammacovid[5],vecsigmacovid[5]))

ysars=log(ddgpd(0:100,vecgammasars[3],vecsigmasars[3]))

ycovid6=log(ddgpd(0:100,vecgammacovid[1],vecsigmacovid[1]))

ysars6=log(ddgpd(0:100,vecgammasars[1],vecsigmasars[1]))

data <- data.frame(
  Disease = c( rep("SARS-CoV ", 101), rep("SARS-CoV-2", 101)),
  Tail_Index = c(ysars,ycovid),
  Tail_Index6 = c(ysars6,ycovid6),
  Threshold = rep((0:100),2)
)

p=ggplot(data,aes(x=Threshold, y=Tail_Index,col=Disease)) +
  scale_color_manual(values = c("red", "blue")) +
  geom_point(aes(x=Threshold, y=Tail_Index,col=Disease),size=0.1) + geom_line() +
  geom_point(aes(x=Threshold, y=Tail_Index6,col=Disease),size=0.1) + geom_line() +
  ylim(c(-10,0)) + xlim(c(0,101)) +
  geom_line() + 
  #geom_abline(intercept = 0.25, slope = 0,color="red",linetype="dotted") +
  #geom_abline(intercept = 0.8, slope = 0,color="blue",linetype="dotted") +
  #geom_abline(intercept = 0.6, slope = 0,color="black",linetype="dotted") + 
  #annotate(geom="text",x=22, y=0.8, label=expression(paste(xi,'=0.8')),color="blue") +
  #annotate(geom="text",x=22, y=0.6, label=expression(paste(xi,'=0.6')),color="black") +
  #annotate(geom="text",x=22, y=0.25, label=expression(paste(xi,'=0.25')),color="red") +
  ylab("log(Probability)") + xlab("x")


print(p+labs(title = "(f) D-GPD log-mass functions")#+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      #                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      ,vp=define_region(3,2))