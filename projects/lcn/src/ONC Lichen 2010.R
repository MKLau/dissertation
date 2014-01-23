#ONC Lichen Data: collected 04May10 by L.J. Lamit and M.K. Lau. 10 cm^2 quadrats on the North and South of each tree placed at 0-10, 45-55, 80-90 cm. Lichen species were observed by L.J. Lamit, data was recorded by M.K. Lau. Lichen identifications were based on identifications of sample specimens by Rikke Naesborg.


#Co-Occurrence Analysis
source('/Users/aeolus/Documents/R_Docs/Scripts/se.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/CU function.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/simIS.R')
source('/Users/aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
L2010=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC Lichens/ONCLichenCooc2010.csv')
setwd('/Users/Aeolus/Documents/Active_Projects/ONC Lichens/')
L.data=read.csv('LCO_03May2011.csv') #lichen co-occur

attach(L2010)
colnames(L2010)

barplot((table(Geno)/50)[order(table(Geno)/50,decreasing=TRUE)])
sum(table(Geno)/50)
summary(Tree)

barplot(apply(L2010[,5:ncol(L2010)],2,sum),las=2)

L.list=list()

Lichen=L2010[,-1:-4]

for (i in 1:nlevels(Tree)){
	L.list[[i]]=Lichen[Tree==levels(Tree)[i],]
	}

names(L.list)=levels(Tree)

geno=character()

for (i in 1:length(L.list)){
	geno[i]=as.character(Geno[Tree==names(L.list)[i]][1])
	}



#CU for every tree
L.cu=list()
for (i in 1:length(L.list)){
	L.cu[[i]]=cu.pairs(t(L.list[[i]]))
	}
names(L.cu)=names(L.list)

#Observed C Score for every tree
cscore.=numeric()
for (i in 1:length(names(L.cu))){
	cscore.[i]=c.score(t(L.list[[i]]))
	}
	
#SES for every tree
library(vegan)

#sequential swap
burn=500
nits=10000

c.obs=0
ses=0
pval.u=0
pval.l=0


for (i in 1:length(L.list)){
x=L.list[[i]]
x=t(x)
L.cscore=cscore.iswap(x,nits=nits,burn=burn)
c.obs[i]=c.score(x)
ses[i]=(c.score(x)-mean(L.cscore$c.sim))/sd(L.cscore$c.sim)
pval.u[i]=length(L.cscore$c.sim[L.cscore$c.sim>=c.obs[i]])/length(L.cscore$c.sim)
pval.l[i]=length(L.cscore$c.sim[L.cscore$c.sim<=c.obs[i]])/length(L.cscore$c.sim)
if (is.na(ses[i])){ses[i]=0}else{}
	}

tree=names(L.list)
out=data.frame(tree,geno,c.obs,ses,pval.u,pval.l)
#write.csv(out,'LCO_03May2011.csv',row.names=FALSE) #lichen co-occur

library(gplots)
attach(L.data)
ses.mu=tapply(ses,geno,mean)
ses.se=tapply(ses,geno,se)

ses.se=ses.se[order(ses.mu)]
ses.mu=ses.mu[order(ses.mu)]

barplot2(ses.mu,plot.ci=TRUE,ci.u=ses.mu+ses.se,ci.l=ses.mu-ses.se)
summary(aov(ses~factor(geno)))

total.A=unlist(lapply(L.list,sum))
plot(total.A~factor(geno))
summary(aov(total.A~factor(geno)))

plot(ses~total.A)
summary(lm(ses~total.A))
summary(lm(ses~total.A*factor(geno)))


L.cu
L.heat=L.cu[[1]]*0

for (i in 1:length(L.cu)){
	L.heat=L.heat+L.cu[[i]]
	}

image(L.heat)

#Network graphs
quartz('Lichen CU Pairs Networks',11,5)
par(mfrow=c(3,7))
for (i in 1:length(L.cu)){gplot(L.cu[[i]],gmode='graph',main=geno[i])}

gplot(L.heat,gmode='graph',main='Sum of CU Pairs')

#scale nodes by sum of CU pairs
v.cex=apply(L.heat,2,sum)/sum(apply(L.heat,2,sum))*15
gplot(L.heat,gmode='graph',main='Sum of CU Pairs',vertex.cex=v.cex,vertex.sides=100,displaylabels=TRUE,vertex.col='gray',edge.col='gray',mode='circle')

#scale nodes by relative species totals 
rst=apply(L2010[,-1:-4],2,sum)/sum(apply(L2010[,-1:-4],2,sum)) #relative species totals
v.cex=rst*10
gplot(L.heat,gmode='graph',main='Sum of CU Pairs',vertex.cex=v.cex,vertex.sides=100,displaylabels=TRUE,vertex.col='gray',edge.col='gray',mode='circle')

#scale nodes by centrality
detach(package:bipartite)
detach(package:tnet)
detach(package:igraph)
v.cex=degree(L.heat)/sum(degree(L.heat))*50
gplot(L.heat,gmode='graph',main='Sum of CU Pairs',vertex.cex=v.cex,vertex.sides=100,displaylabels=TRUE,vertex.col='gray',edge.col='gray',mode='circle')

detach(L.data)

#selecting more trees
tree.list=read.csv('/Users/Aeolus/Documents/Active_Projects/ONC Lichens/LichenTreeList.csv')
geno.list=tree.list$geno
tree.list=tree.list$Tree

geno.n=table(geno.list)/4
geno.names=names(geno.n)[geno.n>=4]
sample(geno.names,1)

unique(tree.list[geno.list=='10'])
