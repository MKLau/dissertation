#Lichen Co-occurrence study
#Sampled the presence of all species of known ONC lichen in the N45-55 lichen quadrats.


###TO DO###
##1. Check that the trees are matched with the appropriate genotypes

##Check for an effect of year?
##Check for an effect of genotype
##Calculate the SES values using the fixed-equiprobable algorithm and the C-score


#Correlation networks
library(vegan)

library(gdata)
library(gplots)
library(sna)
library(lme4)
library(RLRsim)
source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')

##New functions

##Coerce tree notation
as.mytree <- function(x){

x <- as.character(x) #change from factor to character
x <- gsub('-','.',x) #change - to .

                                        #remove leading zeros in the tree number
for (i in (1:length(x))){
  y <- unlist(strsplit(x[i],split='\\.'))
  z <- unlist(strsplit(y[2],split=''))
  if (z[1] == '0'){
    x[i] <- paste(y[1],z[2],sep='.')
  }
}

return(x)
}

##Load data
lco=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')

sep.quad <- TRUE #separate quadrats (N45 and N80)
#lco=lco[lco$Quadrat=='n45.55',]#use only the N45-55 quadrats

summary(lco)
LCO <- lco
if (sep.quad == TRUE){
  LCO$Tree=paste(LCO$Tree,LCO$Quadrat,sep='_')
  tree.quad <- LCO$Tree
} #separate quadrats
names(LCO)

lco.nets=list()

for (i in 1:length(unique(LCO$Tree))){
  lco.nets[[i]]=kendall.pairs(LCO[LCO$Tree==unique(LCO$Tree)[i],7:15],p.adj=FALSE,adj.method='fdr',alpha=0.1) #correlation nets
  ##lco.nets[[i]]=araujo.pairs(LCO[LCO$Tree==unique(LCO$Tree)[i],7:15]) #araujo method nets
}

names(lco.nets) <- unique(LCO$Tree) #names trees using tree number and quadrat

if (sep.quad == TRUE){
lco=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')
LCO <- lco
}

##NA values <- 0
for (i in 1:length(lco.nets)){
  lco.nets[[i]][is.na(lco.nets[[i]])]=0
}

##Re-name each network using the tree name
if (sep.quad == TRUE){
tq.names <- unlist(strsplit(as.character(unique(names(lco.nets))),split='_'))
names(lco.nets) <- tq.names[which(1:length(tq.names) %% 2 == 1)]
}

#genotype factor for each network
geno.net=character()
for (i in 1:length(lco.nets)){
  geno.net[i]=as.character(LCO$Geno[LCO$Tree==names(lco.nets)[i]][1])
}

#quad factor for each network
quad.net=character()
if (sep.quad == TRUE){
  quad.net <- unique(tree.quad)
}else{
for (i in 1:length(lco.nets)){
  quad.net[i]=as.character(LCO$Quadrat[LCO$Tree==names(lco.nets)[i]][1])
}
}

tree.net=character()
for (i in 1:length(lco.nets)){
  tree.net[i]=as.character(lco$Tree[LCO$Tree==names(lco.nets)[i]][1])
}

net.deg=0
for (i in 1:length(lco.nets)){
  net.deg[i]=length(lco.nets[[i]][upper.tri(lco.nets[[i]])][lco.nets[[i]][upper.tri(lco.nets[[i]])]!=0])
}

net.c=0
for (i in 1:length(lco.nets)){
  net.c[i]=centralization(lco.nets[[i]],'degree')
}



##Extract the roughness data for the trees
##For now, use the average roughness from the two quadrats
rough <- read.csv('data/roughness_ONC.csv')
colnames(rough)

attach(rough)

bark.net <- numeric(length=length(tree.net))

Tree <- as.mytree(Tree.O)

for (i in 1:length(tree.net)){
  bark.net[i] <- Roughness.O[Tree == tree.net[i]]
}

bark.net
detach(rough)


plot(geno.net,bark.net)

d.net <- as.dist(netDist(lco.nets))
adonis(d.net~geno.net)
adonis(d.net~bark.net)
adonis(d.net~geno.net*bark.net,perumations=9999)


summary(aov(bark.net~geno.net))


#Accross the entire ONC
lco.sum=lco.nets[[1]]
for (i in 2:length(lco.nets)){
	lco.sum=lco.sum+lco.nets[[i]]
	}
lco.mu=lco.sum/length(lco.nets) #Average tau for all species pairs
gplot(abs(lco.mu),displaylabels=TRUE,gmode='graph')

##For presentation
net.list <- list(lco.mu)

quartz('',15,15)
par(mfrow=c(1,1),mar=c(2.4,1.3,1.3,1.3),oma=c(0.1,0.1,0.1,0.1),bg='transparent',col.main='white',cex=2,mar=c(2,1,1,1))

i <- 1
v.col=apply(net.list[[i]],2,sum); v.col[v.col != 0] = 'lightblue' #color the present species
v.col[v.col == 0] <- 'black' #color empty species
gplot(abs(net.list[[i]]),gmode='graph',vertex.cex=3,vertex.sides=100,vertex.col=v.col,edge.lwd=0.35,edge.col=gray(0.9)[1],vertex.border='grey',mode='circle',displaylabels=TRUE,cex=2,main=names(net.list)[i],label.col='white') #without titles



#Genotype mean networks
geno.munet=tapply(lco.nets,geno.net,edge.mean)

quartz('',10,7)
par(mfrow=c(4,4),mar=c(5, 4, 4, 2)-1)
for (i in (1:length(geno.munet))){
	gplot(geno.munet[[i]],gmode='graph',mode='circle',vertex.sides=1000,vertex.border='black',vertex.col='black',edge.col='grey',main=names(geno.munet)[i],edge.lwd=((geno.munet[[i]]+10)*2))
	}

unique(G)
names(geno.munet)



###For presentaiton
net.list <- geno.munet

quartz('',15,15)
par(mfrow=c(3,5),mar=c(2.4,1.3,1.3,1.3),oma=c(0.1,0.1,0.1,0.1),bg='transparent',col.main='white',cex=2,mar=c(2,1,1,1))

for (i in (1:length(net.list))){
          v.col=apply(net.list[[i]],2,sum); v.col[v.col != 0] = 'lightblue' #color the present species
                              v.col[v.col == 0] <- 'black' #color empty species
                              gplot(abs(net.list[[i]]),gmode='graph',vertex.cex=3,vertex.sides=100,vertex.col=v.col,edge.lwd=0.35,edge.col=gray(0.9)[1],vertex.border='grey',mode='circle',displaylabels=FALSE,cex=2,main=names(net.list)[i]) #without titles
        }



#number of edges
n.edge = numeric()

for (i in 1:length(lco.nets)){
  n.edge[i] = length(lco.nets[[i]][lco.nets[[i]]!=0])
}
  
#abundance of X. gal
XG=numeric()

for (i in 1:length(unique(LCO$Tree))){
  XG[i]=sum(LCO[LCO$Tree==unique(LCO$Tree)[i],7])
}


#percent cover
pcov = numeric()

for (i in 1:length(unique(LCO$Tree))){
  x=apply(LCO[LCO$Tree==unique(LCO$Tree)[i],7:15],2,sum)
  x[x!=0] = 1
  pcov[i]= sum(x) / length(x)
}

pcov.mu=tapply(pcov,geno.net,mean)
pcov.se=tapply(pcov,geno.net,function(x) sd(x)/sqrt(length(x)))
barplot2(pcov.mu,plot.ci=TRUE,ci.u=pcov.mu+pcov.se,ci.l=pcov.mu-pcov.se,las=2)

plot(pcov,n.edge)
abline(lm(n.edge~pcov))
summary(lm(n.edge~pcov))



geno.net=factor(geno.net)
quad.net=factor(quad.net)
tree.net=factor(tree.net)
sqrt.netc=sqrt(net.c)
sqrt.netd=sqrt(net.deg)
quad.net.=as.numeric(quad.net)



summary(aov(net.c~geno.net))
summary(aov(net.deg~geno.net))

summary(aov(sqrt(asin(net.c))~geno.net))
summary(aov(sqrt(net.deg)~geno.net))


exactRLRT(lmer(net.c~(1|geno.net)))
exactRLRT(lmer(net.deg~(1|geno.net)))

#summary(aov(net.c~geno.net*quad.net+Error(tree.net/(quad.net))))
#summary(aov(sqrt.netc~geno.net*quad.net+Error(tree.net/(quad.net))))
#
#summary(aov(net.deg~geno.net*quad.net+Error(tree.net/(quad.net))))
#summary(aov(sqrt.netd~geno.net*quad.net+Error(tree.net/(quad.net))))

#REML
int.net=rep(1,length(sqrt.netc))
full.net=lmer(sqrt.netc~(1|geno.net)+(1|tree.net)+(1|quad.net)+(1|int.net))
gfit.net=lmer(sqrt.netc~(1|geno.net)+(1|int.net))
qtfit.net=lmer(sqrt.netc~(1|tree.net)+(1|quad.net)+(1|int.net))
qfit.net=lmer(sqrt.netc~(1|quad.net)+(1|int.net))
intfit.net=lmer(sqrt.netc~(1|int.net))
anova(full.net,gfit.net,qtfit.net,qfit.net,intfit.net)

int.net=rep(1,length(sqrt.netd))
full.net=lmer(sqrt.netd~(1|geno.net)+(1|tree.net)+(1|quad.net)+(1|int.net))
gfit.net=lmer(sqrt.netd~(1|geno.net)+(1|int.net))
qtfit.net=lmer(sqrt.netd~(1|tree.net)+(1|quad.net)+(1|int.net))
qfit.net=lmer(sqrt.netd~(1|quad.net)+(1|int.net))
intfit.net=lmer(sqrt.netd~(1|int.net))
anova(full.net,gfit.net,qtfit.net,qfit.net,intfit.net)


#The whole garden.
LCO=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/se.R')

x=LCO[,7:15]
colnames(x)
c.sim=cscore.iswap(t(x),nits=10000,burn=500,prow=TRUE,pcol=FALSE)$c.sim


#Using both N45-55 and N80-90 quadrats.
LCO=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/se.R')
summary(LCO)
LCO[LCO$Tree==unique(LCO$Tree)[i],7:14]
barplot(table(LCO$Tree),las=2)

nits=5000
burn=500
tree=character()
geno=character()
year=numeric()
c.obs=numeric()
ses=numeric()
pval.l=numeric()
pval.u=numeric()

date()
play(sin(1:10000/10))
for (i in 1:length(unique(LCO$Tree))){

	tree[i]=as.character(LCO$Tree[LCO$Tree==unique(LCO$Tree)[i]][1])
	geno[i]=as.character(LCO$Geno[LCO$Tree==unique(LCO$Tree)[i]][1])
	year[i]=LCO$Year[LCO$Tree==unique(LCO$Tree)[i]][1]

	x=LCO[LCO$Tree==unique(LCO$Tree)[i],7:15]
	
	c.obs[i]=c.score(t(x))
	c.sim=cscore.iswap(t(x),nits=nits,burn=burn,prow=TRUE,pcol=FALSE)$c.sim
	
	if ((c.obs[i]-mean(c.sim))==0){ses[i]=0}else {ses[i]=(c.obs[i]-mean(c.sim))/sd(c.sim)}
	
	pval.l[i]=length(c.sim[c.sim<=c.obs[i]])/length(c.sim)
	pval.u[i]=length(c.sim[c.sim>=c.obs[i]])/length(c.sim)
	print((i/length(unique(LCO$Tree))))
	}
date()
play(sin(1:10000/10))
out=data.frame(tree,geno,year,c.obs,ses,pval.l,pval.u)
setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
write.csv(out,paste('LCO_results',date(),'.csv',sep=''))

plot(out$ses~out$geno,las=2)
library(nlme)
library(car)

fit.geno=aov(ses~geno,data=out)
summary(fit.geno)
plot(residuals(fit.geno)~fitted(fit.geno))

leveneTest(ses~geno,data=out)

plot(ses~as.numeric(geno),data=out,xaxt='n',xlab='Genotype',ylab='SES')
axis(1,at=unique(as.numeric(out$geno)),labels=levels(out$geno),las=2)

#Percent significant segregated and aggregated
pc.SG=ses #percent segregated/aggregated
pc.SG[pval.l>0.05]=0
pc.SG[pc.SG!=0]=1

pc=function(x){
	x[x!=0]=1
	as.numeric((length(x)-table(x)[1])/length(x))
	}

barplot(tapply(pc.SG,geno,pc),las=2)


#Separate N45-55 and N80-90

LCO=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/se.R')
summary(LCO)
barplot(table(LCO$Tree),las=2)

nits=10000
burn=500
tree=character()
geno=character()
year=numeric()
c.obs=numeric()
ses=numeric()
pval.l=numeric()
pval.u=numeric()

LCO$Tree=paste(LCO$Tree,LCO$Quadrat,sep='_')

date()
play(sin(1:10000/10))
for (i in 1:length(unique(LCO$Tree))){

	tree[i]=as.character(LCO$Tree[LCO$Tree==unique(LCO$Tree)[i]][1])
	geno[i]=as.character(LCO$Geno[LCO$Tree==unique(LCO$Tree)[i]][1])
	year[i]=LCO$Year[LCO$Tree==unique(LCO$Tree)[i]][1]
	
	x=LCO[LCO$Tree==unique(LCO$Tree)[i],7:15]
	
	c.obs[i]=c.score(t(x))
	c.sim=cscore.iswap(t(x),nits=nits,burn=burn,prow=TRUE,pcol=FALSE)$c.sim
	
	if ((c.obs[i]-mean(c.sim))==0){ses[i]=0}else {ses[i]=(c.obs[i]-mean(c.sim))/sd(c.sim)}
	
	pval.l[i]=length(c.sim[c.sim<=c.obs[i]])/length(c.sim)
	pval.u[i]=length(c.sim[c.sim>=c.obs[i]])/length(c.sim)
	print((i/length(unique(LCO$Tree))))
	}
date()
play(sin(1:10000/10))
out=data.frame(tree,geno,year,c.obs,ses,pval.l,pval.u)
setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
write.csv(out,paste('LCO_results',date(),'.csv',sep=''),row.names=FALSE)

lco.data=read.csv("LCO_resultsSun May 15 18:05:51 2011.csv")
tree=as.character(lco.data$tree)

quadrat=as.character(sapply(tree,function(x) strsplit(x,split='_')[[1]][2]))
tree=as.character(sapply(tree,function(x) strsplit(x,split='_')[[1]][1]))
ses=lco.data$ses
geno=lco.data$geno

quad.=as.numeric(factor(quadrat))
tree.=(factor(tree))

summary(aov(sqrt(ses+10)~geno*quad.+Error(tree./(quad.))))

se=function(x){sd(x)/sqrt(length(x))}
mu=tapply(lco.data$ses,lco.data$geno,mean)
se=tapply(lco.data$ses,lco.data$geno,se)
barplot2(mu,las=2,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='SES')

quadrat=sapply(lco.data$tree,)

aov(ses~geno+)

par(mfrow=c(1,2))
plot(ses[quadrat=='n45.55']~geno[quadrat=='n45.55'],las=2)
plot(ses[quadrat=='n80.90']~geno[quadrat=='n80.90'],las=2)

tree.n45=tree[quadrat=='n45.55']
ses.n45=ses[quadrat=='n45.55'][order(tree.n45)]
geno.n45=geno[quadrat=='n45.55'][order(tree.n45)]
tree[quadrat=='n45.55'][order(tree.n45)]

tree.80=tree[quadrat=='n80.90']
ses.n80=ses[quadrat=='n80.90'][order(tree.80)]
geno.n80=geno[quadrat=='n80.90'][order(tree.80)]
tree.80[order(tree.80)]
tree.80[order(tree.80)]==tree[quadrat=='n45.55'][order(tree.n45)]

tree.=tree.80[order(tree.80)]
geno.=geno.n80
ses.=ses.n45-ses.n80

summary(aov(ses.~geno.))

#elevation
elev=as.character(geno.)
elev[elev=='WC5'|elev=='RL6'|elev=='T6']=2
elev[elev!=2]=1
elev=factor(elev)

summary(aov(ses.~elev/geno.))

#REML of ses
setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
lco.data=read.csv("LCO_resultsSun May 15 18:05:51 2011.csv")
names(lco.data)
intercept=rep(1,length(lco.data$ses))
tree.=strsplit(as.character(lco.data$tree),split='_')
tree=character()
for (i in (1:length(tree.))){
	tree[i]=tree.[[i]][1]
	}
quad=character()
for (i in (1:length(tree.))){
	quad[i]=tree.[[i]][2]
	}
geno=lco.data$geno

ses.fit=lmer(ses~(1|geno)+(1|tree)+(1|quad)+(1|intercept))
ses.g.fit=lmer(ses~(1|geno)+(1|intercept))
ses.qt.fit=lmer(ses~(1|tree)+(1|quad)+(1|intercept))
ses.q.fit=lmer(ses~(1|quad)+(1|intercept))
ses.int.fit=lmer(ses~(1|intercept))
anova(ses.fit,ses.g.fit,ses.qt.fit,ses.q.fit,ses.int.fit)

ses.dif=tapply(ses,tree,function(x) x[1]-x[2])
geno.dif.=tapply(geno,tree,function(x) x)
geno.dif=character()
for (i in (1:length(geno.dif.))){
	geno.dif[i]=as.character(geno.dif.[[i]][1])
	}
geno.dif=factor(geno.dif)	
int.dif=rep(1,length(ses.dif))

sesdif.fit=lmer(ses.dif~(1|geno.dif)+(1|int.dif))
intdif.fit=lmer(ses.dif~(1|int.dif))
anova(sesdif.fit,intdif.fit)

log.sesdif=log(ses.dif+10)
logsesdif.fit=lmer(log.sesdif~(1|geno.dif)+(1|int.dif))
logintdif.fit=lmer(log.sesdif~(1|int.dif))
anova(logsesdif.fit,logintdif.fit)


#Using only the N45-55 quadrats. This was collected prior to the N80-90 quadrats.

source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/se.R')
LCO=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_May2011FINAL.csv')
summary(LCO)

nits=10000
burn=500
tree=character()
geno=character()
year=numeric()
c.obs=numeric()
ses=numeric()
pval.l=numeric()
pval.u=numeric()


for (i in 1:length(unique(LCO$Tree))){
	date()
	tree[i]=as.character(LCO$Tree[LCO$Tree==unique(LCO$Tree)[i]][1])
	geno[i]=as.character(LCO$Geno[LCO$Tree==unique(LCO$Tree)[i]][1])
	year[i]=LCO$Year[LCO$Tree==unique(LCO$Tree)[i]][1]
	
	x=LCO[LCO$Tree==unique(LCO$Tree)[i],6:14]
	
	c.obs[i]=c.score(t(x))
	c.sim=cscore.iswap(t(x),nits=nits,burn=0)$c.sim
	
	if ((c.obs[i]-mean(c.sim))==0){ses[i]=0}else {ses[i]=(c.obs[i]-mean(c.sim))/sd(c.sim)}
	
	pval.l[i]=length(c.sim[c.sim<=c.obs[i]])/length(c.sim)
	pval.u[i]=length(c.sim[c.sim>=c.obs[i]])/length(c.sim)
	print((i/length(unique(LCO$Tree))))
	}
date()
out=data.frame(tree,geno,year,c.obs,ses,pval.l,pval.u)
setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
write.csv(out,paste('LCO_results',date(),'.csv',sep=''))

setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
lco.data=read.csv('LCO_resultsSun May  8 20:06:47 2011.csv')[,-1]

geno=factor(geno)
year=factor(year)



plot(ses~geno)
library(gplots)
barplot2(tapply(ses,geno,mean),plot.ci=TRUE,ci.u=tapply(ses,geno,mean)+tapply(ses,geno,se),ci.l=tapply(ses,geno,mean)-tapply(ses,geno,se))




#REML
library(lme4)
library(RLRsim)

lco.reml<-lmer(ses~(1|geno),REML=TRUE,data=lco.data)
summary(lco.reml)
exactRLRT(lco.reml)

#Other calculation of p-value
intercept=rep(1,length(lco.data$ses)) #intercept 
intercept.lmer= lmer(ses~(1|intercept),REML=TRUE,data=lco.data) #intercept only model
lco.lmer<-lmer(ses~(1|geno)+(1|intercept),REML=TRUE,data=lco.data) #full model
anova(intercept.lmer,lco.lmer)
fligner.test(ses~geno,data=lco.data)

plot(ses~geno,data=lco.data)
plot(lco.data$ses~as.numeric(geno),ylab='SES',xlab='Genotype',xaxt='n')
axis(1,labels=levels(geno),at=1:length(levels(geno)),las=2)
plot(tapply(lco.data$ses,lco.data$geno,sd)~factor(levels(lco.data$geno)),las=2)


#kick out H10
ses.noH10=lco.data$ses[lco.data$geno!='H10']
geno.noH10=as.character(lco.data$geno[lco.data$geno!='H10'])
geno.noH10=factor(geno.noH10)
noH10.reml<-lmer(ses.noH10~(1|geno.noH10),REML=TRUE)
exactRLRT(noH10.reml)

#kick out trees N2.3 and N2.4, which were an outliers
ses.noN2.34=lco.data$ses[lco.data$tree!='N2.3'&lco.data$tree!='N2.3']
geno.noN2.34=lco.data$geno[lco.data$tree!='N2.3'&lco.data$tree!='N2.3']
noN2.34.lmer=lmer(ses.noN2.34~(1|geno.noN2.34),REML=TRUE)
exactRLRT(noN2.34.lmer)
summary(aov(ses.noN2.34~geno.noN2.34))
table(lco.data$geno)

#kick out all highly variable genotypes
ses.hv=lco.data$ses[lco.data$geno!='11'&lco.data$geno!='1017'&lco.data$geno!='10']
geno.hv=lco.data$geno[lco.data$geno!='11'&lco.data$geno!='1017'&lco.data$geno!='10']
hv.lmer=lmer(ses.hv~1|geno.hv)
exactRLRT(hv.lmer)

#high vs low
elev=as.character(lco.data$geno)
elev[elev=='WC5'|elev=='RL6'|elev=='T6']=2
elev[elev!=2]=1
elev=as.numeric(elev)
full.elev.lmer=lmer(ses~(1|geno)+(1|elev)+(1|intercept),REML=TRUE,data=lco.data)
elev.elev.lmer=lmer(ses~(1|elev)+(1|intercept),REML=TRUE,data=lco.data)
geno.elev.lmer=lmer(ses~(1|geno)+(1|intercept),REML=TRUE,data=lco.data)
I.elev.lmer=lmer(ses~(1|intercept),REML=TRUE,data=lco.data)
anova(I.elev.lmer,geno.elev.lmer,full.elev.lmer)

ses2010=lco.data$ses[lco.data$year=='2010']
geno2010=lco.data$geno[lco.data$year=='2010']

summary(aov(ses2010~geno2010))

ses2011=lco.data$ses[lco.data$year=='2011']
geno2011=lco.data$geno[lco.data$year=='2011']

summary(aov(ses2011~geno2011))

#garden coordinates
LatLon=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/NONCgpsData.csv')[,1:3]

lco.latlon=array(NA,c(nrow(lco.data),2))
colnames(lco.latlon)=colnames(LatLon)[2:3]
for (i in seq(along=lco.data$tree)){
	lco.latlon[i,]=as.numeric(LatLon[LatLon[,1]==as.character(lco.data$tree[i]),2:3])
	}
lco.data$tree[c(3:4,9,23)]
lco.latlon[c(3:4,9,23),1]=c(-111.9995,-111.9995,-111.9996,-111.9995)
lco.latlon[c(3:4,9,23),2]=c(41.24689,41.24704,41.24842,41.24807)
lco.latlon=data.frame(lco.latlon)
SES=lco.data$ses
GENO=lco.data$geno
lco=data.frame(SES,GENO,lco.latlon)
intercept=rep(1,length=nrow(lco))
lcoLL.lmer=lmer(SES~(1|GENO)+(1|Latitude)+(1|Longitude)+(1|intercept),REML=TRUE,data=lco)
lcoI.lmer=lmer(SES~(1|intercept),REML=TRUE,data=lco)
anova(lcoI.lmer,lcoLL.lmer)

summary(aov(SES~Latitude+Longitude+GENO,data=lco))

#glm

exactRLRT(glmer(ses~1|geno,data=lco.data,family=))


#Restrict to screen genotypes
#996, 1005, 1007, 1008, 1012, 1017

#SES for Screen Sites only
setwd('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/LCO_results')
lco.data=read.csv("LCO_resultsSun May 15 18:05:51 2011.csv")
names(lco.data)

intercept=rep(1,length(lco.data$ses))
tree.=strsplit(as.character(lco.data$tree),split='_')
tree=character()
for (i in (1:length(tree.))){
	tree[i]=tree.[[i]][1]
	}
quad=character()
for (i in (1:length(tree.))){
	quad[i]=tree.[[i]][2]
	}

lco.data=lco.data[quad!='n80.90',]#remove N80-90

intercept=rep(1,length(lco.data$ses))
tree.=strsplit(as.character(lco.data$tree),split='_')
tree=character()
for (i in (1:length(tree.))){
	tree[i]=tree.[[i]][1]
	}
quad=character()
for (i in (1:length(tree.))){
	quad[i]=tree.[[i]][2]
	}

geno=lco.data$geno
ses=lco.data$ses

ses=ses[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
geno=geno[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
geno=factor(as.character(geno))
tree=tree[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
quad=quad[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']

summary(aov(ses~geno))

summary(aov(sqrt(abs(ses))~geno))

plot(((ses))~geno)
plot((abs(ses))~geno)
plot(sqrt(abs(ses))~geno)

barplot2(tapply(ses,geno,mean))

#screen centralization
c.scr=net.c[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
g.scr=factor(as.character(geno.net[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']))
q.scr=quad.net[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
t.scr=tree.net[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']
summary(aov(c.scr~g.scr*q.scr+Error(t.scr/(q.scr))))

mu.c.scr=rbind(tapply(c.scr[q.scr=='n45.55'],g.scr[q.scr=='n45.55'],mean),tapply(c.scr[q.scr=='n80.90'],g.scr[q.scr=='n80.90'],mean))
se.c.scr=rbind(tapply(c.scr[q.scr=='n45.55'],g.scr[q.scr=='n45.55'],se),tapply(c.scr[q.scr=='n80.90'],g.scr[q.scr=='n80.90'],se))

barplot2(mu.c.scr,plot.ci=TRUE,ci.u=mu.c.scr+se.c.scr,ci.l=mu.c.scr-se.c.scr,beside=TRUE,ylab='Centralization')

#screen size
s.scr=net.deg[geno.net=='996'|geno.net=='1005'|geno.net=='1007'|geno.net=='1008'|geno.net=='1012'|geno.net=='1017']

summary(aov(s.scr~g.scr*q.scr+Error(t.scr/(q.scr))))

mu.s.scr=rbind(tapply(s.scr[q.scr=='n45.55'],g.scr[q.scr=='n45.55'],mean),tapply(s.scr[q.scr=='n80.90'],g.scr[q.scr=='n80.90'],mean))
se.s.scr=rbind(tapply(s.scr[q.scr=='n45.55'],g.scr[q.scr=='n45.55'],se),tapply(s.scr[q.scr=='n80.90'],g.scr[q.scr=='n80.90'],se))

barplot2(mu.s.scr,plot.ci=TRUE,ci.u=mu.s.scr+se.s.scr,ci.l=mu.s.scr-se.s.scr,beside=TRUE,ylab='Size')

#screen barplots
quartz("",8,4.5)
par(mfrow=c(1,2))
barplot2(mu.s.scr,plot.ci=TRUE,ci.u=mu.s.scr+se.s.scr,ci.l=mu.s.scr-se.s.scr,beside=TRUE,ylab='Size',col=c('black','gray'),las=2)
legend('topright',legend=c('45-55cm','80-90cm'),fill=c('black','gray'),bty='n')
barplot2(mu.c.scr,plot.ci=TRUE,ci.u=mu.c.scr+se.c.scr,ci.l=mu.c.scr-se.c.scr,beside=TRUE,ylab='Centralization',col=c('black','gray'),las=2)
legend('topright',legend=c('45-55cm','80-90cm'),fill=c('black','gray'),bty='n')

#screen reml
#scree centralization reml
detach(package:RLRsim)
detach(package:nlme)
t.scr=factor(as.character(t.scr))
int.scr=rep(1,length(c.scr))

full.scr=lmer(c.scr~(1|g.scr)+(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(c.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(c.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(c.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(c.scr~1|int.scr,verbose=TRUE)
anova(full.scr,gfit.scr,qfit.scr,qtfit.scr,intfit.scr)

sqrt.c.scr=sqrt(c.scr)
full.scr=lmer(sqrt.c.scr~(1|g.scr)+(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(sqrt.c.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(sqrt.c.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(sqrt.c.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(sqrt.c.scr~1|int.scr,verbose=TRUE)
anova(gfit.scr,full.scr,qfit.scr,qtfit.scr,intfit.scr)

log.c.scr=log(c.scr+0.001)
full.scr=lmer(log.c.scr~(1|g.scr)+(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(log.c.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(log.c.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(log.c.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(log.c.scr~1|int.scr,verbose=TRUE)
anova(gfit.scr,full.scr,qfit.scr,qtfit.scr,intfit.scr)

#screen size reml
full.scr=lmer(s.scr~(1|g.scr)+(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(s.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(s.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(s.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(s.scr~1|int.scr,verbose=TRUE)
anova(gfit.scr,full.scr,qfit.scr,qtfit.scr,intfit.scr)

sqrt.s.scr=sqrt(s.scr)
full.scr=lmer(sqrt.s.scr~(1|t.scr)+(1|q.scr)+(1|g.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(sqrt.s.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(sqrt.s.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(sqrt.s.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(sqrt.s.scr~1|int.scr,verbose=TRUE)
anova(gfit.scr,qfit.scr,qtfit.scr,intfit.scr,full.scr)

log.s.scr=log(s.scr+0.001)
full.scr=lmer(log.s.scr~(1|g.scr)+(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
gfit.scr=lmer(log.s.scr~(1|g.scr)+(1|int.scr),verbose=TRUE)
qfit.scr=lmer(log.s.scr~(1|q.scr)+(1|int.scr),verbose=TRUE)
qtfit.scr=lmer(log.s.scr~(1|t.scr)+(1|q.scr)+(1|int.scr),verbose=TRUE)
intfit.scr=lmer(log.s.scr~1|int.scr,verbose=TRUE)
anova(gfit.scr,full.scr,qfit.scr,qtfit.scr,intfit.scr)

#Difference screen reml
dc.scr=tapply(c.scr,t.scr,function(x) x[1]-x[2])
dg.scr=unlist(tapply(g.scr,t.scr,function(x) return(x)))
dg.scr=factor(dg.scr[(1:length(dg.scr))%%2==1])
dc.int=rep(1,length(dc.scr))

full.scr=lmer(dc.scr~(1|dg.scr)+(1|dc.int),verbose=TRUE)
dc.intfit.scr=lmer(dc.scr~(1|dc.int),verbose=TRUE)
anova(dc.intfit.scr,full.scr)

ds.scr=tapply(s.scr,t.scr,function(x) x[1]-x[2])
dg.scr=unlist(tapply(g.scr,t.scr,function(x) return(x)))
dg.scr=factor(dg.scr[(1:length(dg.scr))%%2==1])
dc.int=rep(1,length(ds.scr))

full.scr=lmer(ds.scr~(1|dg.scr)+(1|dc.int),verbose=TRUE)
dc.intfit.scr=lmer(ds.scr~(1|dc.int),verbose=TRUE)
anova(dc.intfit.scr,full.scr)

#Graphs
#Barplot
library(gplots)

#Average network structure for upper and lower
munet.u=lco.nets[[1]]*0
for (i in ((1:length(quad.net))[quad.net=='n80.90'])){
	munet.u=munet.u+lco.nets[[i]]
	}
munet.u=munet.u/length(((1:length(quad.net))[quad.net=='n80.90'])) #Average tau for all species pairs
gplot(abs(munet.u),displaylabels=TRUE,gmode='graph')

munet.l=lco.nets[[1]]*0
for (i in ((1:length(quad.net))[quad.net=='n45.55'])){
	munet.l=munet.l+lco.nets[[i]]
	}
munet.l=munet.l/length(((1:length(quad.net))[quad.net=='n45.55'])) #Average tau for all species pairs
gplot(abs(munet.l),displaylabels=TRUE,gmode='graph')

par(mfrow=c(1,2))
gplot(abs(munet.l),displaylabels=TRUE,gmode='graph',vertex.col='black',edge.col='gray',xlab='45-55cm',mode='circle')
gplot(abs(munet.u),displaylabels=TRUE,gmode='graph',vertex.col='black',edge.col='gray',xlab='80-90cm',mode='circle')

#Upper versus lower
barplot(tapply(net.deg,quad.net,mean),ylab='Size')
barplot(tapply(net.c,quad.net,mean),ylab='Centralization')

#Size
mu.u=tapply(net.deg[quad.net=='n80.90'],geno.net[quad.net=='n80.90'],mean)
mu.l=tapply(net.deg[quad.net=='n45.55'],geno.net[quad.net=='n45.55'],mean)
mu.d=mu.l-mu.u
mu.u=mu.u[order(abs(mu.d),decreasing=TRUE)]
mu.l=mu.l[order(abs(mu.d),decreasing=TRUE)]

se=function(x){sd(x)/sqrt(length(x))}

se.u=tapply(net.deg[quad.net=='n80.90'],geno.net[quad.net=='n80.90'],se)
se.l=tapply(net.deg[quad.net=='n45.55'],geno.net[quad.net=='n45.55'],se)
se.u=se.u[order(abs(mu.d),decreasing=TRUE)]
se.l=se.l[order(abs(mu.d),decreasing=TRUE)]

par(mfrow=c(1,2))
barplot2(rbind(mu.l,mu.u),beside=TRUE,las=2,col=c('black','gray'),plot.ci=TRUE,ci.u=rbind(mu.l+se.l,mu.u+se.u),ci.l=rbind(mu.l-se.l,mu.u-se.u),ylab='Size',xlab='')
legend('topright',legend=c('45-55cm','80-90cm'),bty='n',fill=c('black','gray'),border='white')


#Centralization
mu.u=tapply(net.c[quad.net=='n80.90'],geno.net[quad.net=='n80.90'],mean)
mu.l=tapply(net.c[quad.net=='n45.55'],geno.net[quad.net=='n45.55'],mean)
#mu.d=mu.l-mu.u
mu.u=mu.u[order(abs(mu.d),decreasing=TRUE)]
mu.l=mu.l[order(abs(mu.d),decreasing=TRUE)]

se=function(x){sd(x)/sqrt(length(x))}

se.u=tapply(net.c[quad.net=='n80.90'],geno.net[quad.net=='n80.90'],se)
se.l=tapply(net.c[quad.net=='n45.55'],geno.net[quad.net=='n45.55'],se)
se.u=se.u[order(abs(mu.d),decreasing=TRUE)]
se.l=se.l[order(abs(mu.d),decreasing=TRUE)]

barplot2(rbind(mu.l,mu.u),beside=TRUE,las=2,col=c('black','gray'),plot.ci=TRUE,ci.u=rbind(mu.l+se.l,mu.u+se.u),ci.l=rbind(mu.l-se.l,mu.u-se.u),ylab='Centralization',xlab='')
legend('topright',legend=c('45-55cm','80-90cm'),bty='n',fill=c('black','gray'),border='white')

source('/Users/Aeolus/Documents/R_Docs/Scripts/Functions/pairs.R')
pairs(lco[,7:15],lower.panel=panel.cor,upper.panel=panel.lm)

#Given a list of adjacency matrices and an indexing vector, plot the indexed graphs
length(lco.nets)
lco.nets
geno.net

igplot=function(x,y,index){
	if (length(y[y==index])>25){print('ERROR: number of plots too great for window.')}else{
	quartz(index,width=9,height=7)
	par(mfrow=c(round(sqrt(length(y[y==index])),0),round(sqrt(length(y[y==index]))+1,0)))
	for (i in (1:length(y))[y==index]){
		gplot(x[[i]],gmode='graph',vertex.sides=1000,vertex.border='black',vertex.col='black',edge.col='grey')
		}
	}
}

for (i in 1:nlevels(geno.net)){
	igplot(lco.nets,geno.net,levels(geno.net)[i])
	}

##hybrid index
## HI <- read.csv('data/ONC_hybrid_index.csv')
## HI <- na.omit(HI)
## HI.index <- as.numeric(HI[,3])
## HI.geno <- as.character(HI[,1])
## HI.geno <- gsub('-','',HI.geno)
## geno.net <- as.character(geno.net)
## lco.hi <- geno.net

## for (i in (1:length(unique(geno.net)))){
## lco.hi[lco.hi == unique(geno.net)[i]] <- HI.index[HI.geno == unique(geno.net)[i]]
## }
## lco.hi <- as.numeric(lco.hi)
##adonis(d.net~lco.hi)
