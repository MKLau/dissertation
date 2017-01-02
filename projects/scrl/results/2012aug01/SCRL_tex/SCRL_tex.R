###################################################
### chunk number 1: 
###################################################

#remove dead trees
lichen[[3]]=lichen[[3]][lichen[[3]][,3]==1,]

com=lichen[[3]][,13:ncol(lichen[[3]])]

attach(lichen[[3]][,-13:-ncol(lichen[[3]])])
summary(lichen[[3]][,-13:-ncol(lichen[[3]])])
env=lichen[[3]][,-13:-ncol(lichen[[3]])]
env=env[,-3]
tree=rep(c(1,0),(nrow(com)/2))
data.frame(lichen[[1]][1:nrow(com),1],tree)


###################################################
### chunk number 2: 
###################################################
summary(com)
summary(env)


###################################################
### chunk number 3: 
###################################################
library(vegan)
library(gplots)
attach(env)


###################################################
### chunk number 4: 
###################################################
par(mfrow=c(1,2))
plot(specaccum(com),lwd=0.5,xlab='Trees Sampled',ylab='Species Number',main='All Trees') #all trees
plot(specaccum(com[tree==1,]),add=FALSE,xlab='Trees Sampled',ylab='Species Number')#resistant trees
abline(h=15,lwd=0.3,lty=2)
plot(specaccum(com[tree==0,]),add=TRUE,col='red')#susceptible trees
abline(h=13,lwd=0.3,lty=2)
legend('bottomright',c('S','R'),lty=1,col=c('black','red'),bty='n')


###################################################
### chunk number 5: 
###################################################
a=apply(com,1,sum)
wilcox.test(a,y=Tree.pairs,alternative='less',paired=TRUE)


###################################################
### chunk number 6: 
###################################################
r=com
r[r!=0]<-1
r=apply(r,1,sum)
wilcox.test(r,y=Tree.pairs,alternative='less',paired=TRUE)



###################################################
### chunk number 7: 
###################################################
h=diversity(com)
wilcox.test(h,y=Tree.pairs,alternative='less',paired=TRUE)


###################################################
### chunk number 8: 
###################################################
se=function(x){sd(x)/sqrt(length(x))}
par(mfrow=c(1,3))
barplot2(tapply(a,tree,mean),xlab='Resistance',ylab='Percent Abundance',plot.ci=TRUE,ci.u=tapply(a,tree,mean)+tapply(a,tree,se),ci.l=tapply(a,tree,mean)-tapply(a,tree,se))
barplot2(tapply(r,tree,mean),xlab='Resistance',ylab='Species Richness',plot.ci=TRUE,ci.u=tapply(r,tree,mean)+tapply(r,tree,se),ci.l=tapply(r,tree,mean)-tapply(r,tree,se))
barplot2(tapply(h,tree,mean),xlab='Resistance',ylab='Shannon Index',plot.ci=TRUE,ci.u=tapply(h,tree,mean)+tapply(h,tree,se),ci.l=tapply(h,tree,mean)-tapply(h,tree,se))


###################################################
### chunk number 9: 
###################################################
source('/Users/artemis/Documents/R_Docs/Scripts/Functions/pairs.R')
pairs(com[apply(com,2,sum)>=5],lower.panel=panel.cor,upper.panel=panel.lmr)


###################################################
### chunk number 10: 
###################################################
d=vegdist(cbind(com,rep(0.01,nrow(com))))
adonis(d~Moth+Tree.pairs/Moth,permutations=99999)


###################################################
### chunk number 11: 
###################################################
library(labdsv)
indsp=indval(com,(tree+1))
summary(indsp)
detach(package:labdsv)


###################################################
### chunk number 12: 
###################################################
nms.x=read.csv('moth_nms.csv')
plot(nms.x,col=(env$Moth+1))
fit=envfit(nms.x,env[,c(-9,-10)],permutations=99999)
plot(fit,col='black')


###################################################
### chunk number 13: 
###################################################
fit


###################################################
### chunk number 14: 
###################################################
pairs(cbind(Moth,Litter..,Big.rocks..,Small.rocks..,Light...average),lower.panel=panel.cor,upper.panel=panel.lmr)
cor.test(Litter..,Moth)
cor.test(Big.rocks..,Moth)
cor.test(Small.rocks..,Moth)
cor.test(Light...average,Moth)


###################################################
### chunk number 15: 
###################################################
par(mfrow=c(1,3))
hist(residuals(lm(a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='untransformed')
log.a=log(a+0.0000001)
hist(residuals(lm(log.a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='log transformed')
sqrt.a=sqrt(a)
hist(residuals(lm(sqrt.a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='sqrt transformed')

summary(lm(a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))
summary(lm(sqrt.a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))


###################################################
### chunk number 16: 
###################################################
par(mfrow=c(1,3))
hist(residuals(lm(r~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='untransformed')
hist(residuals(lm(sqrt(r)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='sqrt transformed')
hist(residuals(lm(r^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='square transformed')

summary(lm(r~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))
summary(lm(r^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))


###################################################
### chunk number 17: 
###################################################
par(mfrow=c(2,2))
hist(residuals(lm(h~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='untransformed')
hist(residuals(lm(h^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='square transformed')
hist(residuals(lm(sqrt(h)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='sqrt transformed')
hist(residuals(lm(log(h+0.0000001)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)),main='log transformed')

summary(lm(h~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))
summary(lm(h^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))


###################################################
### chunk number 18: 
###################################################
adonis(d~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth,permutations=9999)


