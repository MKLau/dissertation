#Rock Lichen data from Sunset Crater collected by Rikke Naesborg and Richard Michalet

dif=function(x){
	out=x[1]
	for (i in 2:length(x)){
		out=out-x[i]
		}
		return(out)
	}


#Distributions of predictors and responses
#Analyses of all factors with abundance, richness, Shannon's Index and community composition

#IS THE QUADRAT DATA PROBLEMATIC???

#Rank abundance plots for the two tree types
#Species accumulation curve within tree types
#Relativized abundance for community analyes
#Handle the spatial and live-dead structure 


setwd('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/lichen & abiotic cover Sunset Crater paired')

lichen=list()
for (i in 1:length(dir())){
	lichen[[i]]=read.csv(dir()[i])
	}
names(lichen)=dir()
names(lichen)

summary(lichen[[3]])

#remove dead trees
lichen[[3]]=lichen[[3]][lichen[[3]][,3]==1,]

com=lichen[[3]][,13:ncol(lichen[[3]])]



attach(lichen[[3]][,-13:-ncol(lichen[[3]])])
summary(lichen[[3]][,-13:-ncol(lichen[[3]])])
env=lichen[[3]][,-13:-ncol(lichen[[3]])]
env=env[,-3]
summary(env)



tree=rep(c(1,0),(nrow(com)/2))
data.frame(lichen[[1]][1:nrow(com),1],tree)

library(vegan)
library(ecodist)

#Species accumulation curves
par(mfrow=c(1,2))
plot(specaccum(com),lwd=0.5,xlab='Trees Sampled',ylab='Species Number',main='All Trees') #all trees
plot(specaccum(com[tree==1,]),add=FALSE,xlab='Trees Sampled',ylab='Species Number')#resistant trees
abline(h=15,lwd=0.3,lty=2)
plot(specaccum(com[tree==0,]),add=TRUE,col='red')#susceptible trees
abline(h=13,lwd=0.3,lty=2)
legend('bottomright',c('S','R'),lty=1,col=c('black','red'),bty='n')

#by individual ***this is kind of wonky because we are working with percentages, not counts
plot(specaccum(round(com,0),method='rarefaction'),lwd=0.5) #all trees
plot(specaccum(round(com[tree==1,],0),method='rarefaction'),add=TRUE)#resistant trees
plot(specaccum(round(com[tree==0,],0),method='rarefaction'),add=TRUE,col='red')#susceptible trees
legend('bottomright',c('S','R','All'),lty=1,col=c('black','red','grey'),bty='n')

#Distributions of Predictors
source('/Users/artemis/Documents/R_Docs/Scripts/Functions/pairs.R')
pairs(env,lower.panel=panel.cor)
pairs(env,lower.panel=panel.hist)

#Abundance, richness and evenness
#The difference is Resistant - Susceptible
a=apply(com,1,sum)
r=com
r[r!=0]<-1
r=apply(r,1,sum)
h=diversity(com)
ap=tapply(a,tp,dif)
rp=tapply(r,tp,dif)
hp=tapply(h,tp,dif)

wilcox.test(a,y=Tree.pairs,alternative='less',paired=TRUE)
wilcox.test(r,y=Tree.pairs,alternative='less',paired=TRUE)
wilcox.test(h,y=Tree.pairs,alternative='less',paired=TRUE)

par(mfrow=c(1,3))
plot(density(a[tree==0]),main='')
lines(density(a[tree==1]),col='red')
legend('topright',c('Sus','Res'),lty=c(1,1),bty='n',col=c(1,2))

plot(density(r[tree==0]),main='')
lines(density(r[tree==1]),col='red')
legend('topright',c('Sus','Res'),lty=c(1,1),bty='n',col=c(1,2))
plot(density(h[tree==1]),col='red',main='')
lines(density(h[tree==0]))
legend('topright',c('Sus','Res'),lty=c(1,1),bty='n',col=c(1,2))

boxplot(a~tree,outline=FALSE,xlab='Resistance',ylab='Abundance')
boxplot(r~tree,outline=FALSE,xlab='Resistance',ylab='Richness')
boxplot(h~tree,outline=FALSE,xlab='Resistance',ylab='Shannon Index')

se=function(x){sd(x)/sqrt(length(x))}

library(gplots)
par(mfrow=c(1,3))
barplot2(tapply(a,tree,mean),xlab='Resistance',ylab='Percent Abundance',plot.ci=TRUE,ci.u=tapply(a,tree,mean)+tapply(a,tree,se),ci.l=tapply(a,tree,mean)-tapply(a,tree,se))
barplot2(tapply(r,tree,mean),xlab='Resistance',ylab='Species Richness',plot.ci=TRUE,ci.u=tapply(r,tree,mean)+tapply(r,tree,se),ci.l=tapply(r,tree,mean)-tapply(r,tree,se))
barplot2(tapply(h,tree,mean),xlab='Resistance',ylab='Shannon Index',plot.ci=TRUE,ci.u=tapply(h,tree,mean)+tapply(h,tree,se),ci.l=tapply(h,tree,mean)-tapply(h,tree,se))

#Multivariate Community Analyses
pairs(com[apply(com,2,sum)>=5],lower.panel=panel.cor,upper.panel=panel.lmr)
#Composition
d=vegdist(cbind(com,rep(0.01,nrow(com))))
adonis(d~Moth+Tree.pairs/Moth,permutations=99999)

#adonis(d~tree,nperm=5000) without pairs

#mrpp(d,tree,strata=Tree.pairs)

#indicator species analysis
library(labdsv)

indsp=indval(com,(tree+1))
summary(indsp)

detach(package:labdsv)

#community ordination
x=com
x=cbind(x,rep(1,nrow(com)))
dx=vegdist(x)
#nms.x=nmds(dx,1,5)
#plot(nms.x$stress)
#nms.x=nmds(dx,2,2,nits=50)
#nms.x=nmds.min(nms.x,2)
#setwd('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/SCRL_tex/')
#write.csv(nms.x,'moth_nms.csv',row.names=FALSE)

#Save ordination
setwd('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/SCRL_tex/')
nms.x=read.csv('moth_nms.csv')
plot(nms.x,col=(env$Moth+1))
fit=envfit(nms.x,env[,c(-9,-10)],permutations=99999)
plot.envfit(fit,col='black')
fit


#Exploring the factors behind the effect of resistance on lichen communities
#Moth and Factors
pairs(cbind(Moth,Litter..,Big.rocks..,Small.rocks..,Light...average),lower.panel=panel.cor,upper.panel=panel.lmr)
cor.test(Litter..,Moth)
cor.test(Big.rocks..,Moth)
cor.test(Small.rocks..,Moth)
cor.test(Light...average,Moth)

#Multiple regression

hist(residuals(lm(a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
log.a=log(a+0.0000001)
hist(residuals(lm(log.a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
sqrt.a=sqrt(a)
hist(residuals(lm(sqrt.a~Moth+Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))

summary(lm(sqrt.a~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))

hist(residuals(lm(r~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
hist(residuals(lm(sqrt(r)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
hist(residuals(lm(r^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))

summary(lm(r^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))

hist(residuals(lm(h~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
hist(residuals(lm(h^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
hist(residuals(lm(sqrt(h)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))
hist(residuals(lm(log(h+0.0000001)~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth)))

summary(lm(h^2~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth))

#Vector Analyses
x=env[,c(-1:-3,-10,-11)]
pairs(x)
dx=vegdist(x)
dx=vegdist(x)
#nms.x=nmds(dx,1,5,nits=10)
#plot(nms.x$stress,type='l')
nms.x=nmds(dx,2,2,nits=100)
nms.x=nmds.min(nms.x,2)
plot(nms.x,col=(env$Moth+1),main='Ordination of Tree via Environmental Parameter')
vec.x=envfit(nms.x,x)
plot(vec.x)
legend('topright',c('S','R'),pch=c(1,1),col=c(1,2))

#PerMANOVA
adonis(d~Litter..+Big.rocks..+Small.rocks..+Light...average+Moth,permutations=9999)


#SEM preliminaries
#correlation structure with abundance, richness and shannon's index
source('/Users/artemis/Documents/R_Docs/Scripts/Functions/pairs.R')
x=cbind(env,a,r,h)
x=x[,c(-1)]
pairs(x,upper.panel=panel.cor,lower.panel=panel.lmr)

x=cbind(env,a)
x=x[,c(-2,-3)]
pairs(x,upper.panel=panel.cor,lower.panel=panel.lmr)

pairs(cbind(a,Litter..,Big.rocks..,Small.rocks..))
plot(Moth,Litter..)

#Write to file for input to AMOS
#SCRL=lichen[[3]]
#SCRL=cbind(SCRL,a,r,h)
#colnames(SCRL)[1:12]=c('pairs','resistance','living','litter','big.rocks','small.rocks','shrubs','grass','branches','light.N','light.S','light.avg')
#summary(SCRL)
#write.csv(SCRL,'SCRL')

#Communities
#com.p=array(NA,c(length(unique(tp)),ncol(com)))
#colnames(com.p)<-colnames(com)
#
#for (i in 1:length(unique(tp))){
#	com.p[i,]=as.matrix(com[tp==unique(tp)[i],][1,]-com[tp==unique(tp)[i],][2,])
#	}
#Correlation Graphs
#library(sna)
#lich.cor=cor(com)
#lich.cor[abs(lich.cor)<0.2]=0
#size=apply(com,2,sum)+1
#gplot(abs(lich.cor),gmode='graph',displaylabels=TRUE,vertex.cex=size/5,vertex.sides=100,vertex.border='red',edge.col='darkgrey')

pairs(com[apply(com,2,sum)>=5])

d=vegdist(cbind(com,rep(0.01,nrow(com))))

adonis(d~tree+Tree.pairs/tree,nperm=5000)
#mrpp(d,tree,strata=Tree.pairs)

#indicator species analysis
library(labdsv)

indsp=indval(com,(tree+1))
summary(indsp)

detach(package:labdsv)

#Ordination
scree=nmds(d,1,3,3)$stress
plot(scree)

ord=nmds(d,3,3,nits=10)
ord=nmds.min(ord,3)
plot(ord,col=(tree+1))
ord.fig=ordiplot(ord,c(1,3),type='n')
points(ord.fig,'sites',col=(tree+1))
ordispider(ord.fig,Tree.pairs)
hist(tree[apply(com,1,sum)==0])
Tree.pairs[apply(com,1,sum)==0]

#Remove the zero sums
my.col=c('green','purple')
ord.rzs=ord[apply(com,1,sum)!=0,]
tree.rzs=tree[apply(com,1,sum)!=0]
Tree.pairs.rzs=Tree.pairs[apply(com,1,sum)!=0]
ord.fig=ordiplot(ord.rzs,c(1,2),type='n')
points(ord.fig,'sites',col=my.col[(tree.rzs+1)],pch=15)
ordispider(ord.fig,Tree.pairs.rzs,col='orange',lwd=0.5)

#
#
#setwd('/Users/artemis/Documents/R_Docs/Rikke Naesborg')
#dir()
#setwd('Rikke Data')
#dir()
#rl.c=read.csv('sp & env combined-Table 1.csv')[,-1:-9]#all observations
#rl.e=read.csv('sp & env combined-Table 1.csv')[,1:9]#all observations
##rl.c=read.csv("Rikke.csv")[,-1:-9] #rock lichen community
##rl.e=read.csv("Rikke.csv")[,1:9] #rock lichen environmental data
#
#library(vegan)
#
#names(rl.e)
#attach(rl.e)
#plot(Moth,Rocks.3cm..)
#
##Explore the abundances
#A=apply(rl.c,2,sum) #total species abundances
#barplot(A[order(A,decreasing=TRUE)],name=colnames(rl.c)[order(A,decreasing=TRUE)],las=2)
#
#Ap=apply(rl.c,1,sum) #total abundance per plot
#
#summary(aov(sqrt(Ap)~factor(Moth)+factor(Live.Dead)+(Litter.)+Rocks.3cm..+Rocks.3cm...1))
#hist(aov(sqrt(Ap)~factor(Moth)+factor(Live.Dead)+(Litter.))$residuals)
#
#plot(Litter.~Rocks.3cm..)
#pairs(rl.e)
#
#hist(aov(sqrt(Ap)~factor(Moth))$residuals)
#barplot(c(mean(Ap[Moth==0]),mean(Ap[Moth==1])),names=c('S','R'),ylab='Total Abundance')
#
##richness and diversity
#quartz()
#R=specnumber(rl.c)
#summary(aov(R~factor(Moth)*factor(Live.Dead)*(Litter.)))
#hist(residuals(aov(R~factor(Moth)+factor(Live.Dead)+(Litter.))))
#
#plot(R~Litter.,col=(as.numeric((Moth))+1))
#abline(lm(R[Moth==0]~Litter.[Moth==0]),col=1)
#abline(lm(R[Moth==1]~Litter.[Moth==1]),col=2)
#
#H=diversity(rl.c)
#summary(aov((H)~factor(Moth)*factor(Live.Dead)*(Litter.)))
#hist(residuals(aov((H)~factor(Moth)+factor(Live.Dead)+(Litter.))))
#
#plot(H~Litter.,col=(as.numeric((Moth))+1))
#abline(lm(H[Moth==0]~Litter.[Moth==0]),col=1)
#abline(lm(H[Moth==1]~Litter.[Moth==1]),col=2)
#
##Removing rare species
#rl.c=rl.c[,A>=5]
#names(rl.c)
#
##check for zero sum rows
#if (any(apply(rl.c,1,sum)==0)){
#dummy=rep(min(rl.c[rl.c!=0]),nrow(rl.c))
#d=vegdist(data.frame(rl.c,dummy))}else{
#d=vegdist(data.frame(rl.c))	
#	}
#
##permanova
#adonis(d~factor(Moth)*factor(Live.Dead)*Litter.)
#
##relativize by species max
##check for zero sum rows
#rl.c.max=decostand(rl.c,2,method='max')
#if (any(apply(rl.c.max,1,sum)==0)){
#dummy=rep(min(rl.c.max[rl.c.max!=0]),nrow(rl.c.max))
#d.max=vegdist(data.frame(rl.c.max,dummy))}else{
#d.max=vegdist(data.frame(rl.c.max))	
#	}
#
##permanova
#adonis(d.max~factor(Moth)*factor(Live.Dead)*Litter.)
#
#library(ecodist)
##n=nmds(vegdist(rl.c),1,5)
##plot(n$stress) 
##three dimensions looks best
#
##ordination
#n=nmds(d,3,3,nits=50)
#plot(n$stress)
#n.=nmds.min(n)
#plot(n.,col=(as.numeric(Moth)+1))
#
#ordiplot(n.,c(1,2),type='none')
#points(n.[,1:2],pch=(as.numeric(Moth)+1))
#ordispider(n.[,1:2],Moth,col=c('violet'))
#plot(envfit(n.~Litter.),add=TRUE)
#
#
##Multivariate Dispersion Test
#b=betadisper(d,Moth,type='centroid')
#permutest.betadisper(b)
#par(mfrow=c(1,2))
#plot(b,main='Dispersion')
#plot(TukeyHSD(b),ylab='')
#
##Indicator Species Analysis
#library(labdsv)
#names(rl.c)
#summary(indval(rl.c,Moth))
#ds=vegdist(t(rl.c))
#plot(hclust(ds))
#
##Cluster diagram of plots with moth overlay
#plot(hclust(d),labels=Moth)
#
#detach(rl.e)
#rm(list=c('d','rl.c','rl.e'))
#
#
##Output to PCORD
##rock.stability=Rocks.3cm../(Rocks.3cm..+Rocks.3cm...1) 
##this index didn't work because sometimes there are only big rocks and no little rocks
##However, in doing this it seems as though using the both factors in the model might be useful
#
#out=cbind(Ap,R,H,rock.stability,rl.e)
#
#write.csv(out,'/Users/artemis/Desktop/rocklichen',row.names=FALSE)
#
####1) Calculate the roack size index (stability index) = A_r>3 / A_rT
####2) Conduct SEM (NOTE! Handle the binary response variables!!!)
#
##Test for the effect of rocks >3cm on moth susceptibility
#rocks.glm=glm(Moth~Rocks.3cm..,family='binomial')
#summary(rocks.glm)
#plot(Moth~Rocks.3cm..)
#points(Rocks.3cm..,fitted(glm(Moth~Rocks.3cm..,family=binomial(link='logit'))),cex=0.5,col='red')
#lines(spline(Rocks.3cm..,fitted(glm(Moth~Rocks.3cm..,family=binomial(link='logit')))),cex=0.5,col='red')
#AIC(glm(Moth~Rocks.3cm..,family=binomial(link='logit')))
#AIC(lm(Moth~Rocks.3cm..))
#points(Rocks.3cm..,fitted(lm(Moth~Rocks.3cm..)),cex=0.5,col='blue')
#lines(Rocks.3cm..,fitted(lm(Moth~Rocks.3cm..)),cex=0.5,col='blue')
#legend('right',legend=c('Linear','Binomial'),bty='n',fill=c('blue','red'),border=c('blue','red'))

#Sweave

setwd('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/lichen & abiotic cover Sunset Crater paired')

lichen=list()
for (i in 1:length(dir())){
	lichen[[i]]=read.csv(dir()[i])
	}
names(lichen)=dir()
names(lichen)

summary(lichen[[3]])

#remove dead trees
lichen[[3]]=lichen[[3]][lichen[[3]][,3]==1,]
my.wd=getwd()
setwd('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/SCRL_tex')
Sweave('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/SCRL_tex/SCRL_tex.Rnw')
Stangle('/Users/artemis/Documents/Active_Projects/Sunset_Crater_Lichens/SCRL_tex/SCRL_tex.Rnw')
setwd(my.wd)
