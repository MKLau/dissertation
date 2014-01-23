#Analysis of arthropod co-occurrence patterns on narrowleaf genotypes. Collected with Art Keith, June 15-17 2011. 

source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')
source('/Users/Aeolus/Documents/R_Docs/Scripts/EcosimR/EcosimR.R')
library(audio)
library(sna)

acn=read.csv('/Users/Aeolus/Documents/Active_Projects/ArthopodCooccurrenceNetworks/ONC_AC_2011.csv')

acn[is.na(acn)]=0

geno=acn$Geno

acn=acn[,-1:-5]
acn[,1]=acn[,11]#remove PB aborts
acn=acn[,-11]#remove PB aborts

acn.net=list()
for (i in (1:length(unique(geno)))){
	acn.net[[i]]=kendall.pairs(acn[geno==unique(geno)[i],])
}

acn.cup=list()
for (i in (1:length(unique(geno)))){
	acn.cup[[i]]=cu.pairs(t(acn[geno==unique(geno)[i],]))
}

names(acn.cup)=unique(geno)

acn.cs=list()
for (i in (1:length(unique(geno)))){
	acn.cs[[i]]=cscore.iswap(t(acn[geno==unique(geno)[i],]),5000,pcol=FALSE)$c.sim
}

play(sin(1:10000/20))

acn.ses=numeric()
for (i in (1:length(unique(geno)))){
	obs=c.score(t(acn[geno==unique(geno)[i],]))
	acn.ses[i]=obs-mean(acn.cs[[i]])/sd(acn.cs[[i]])
}

acn.ses[acn.ses==-Inf]=0

quartz('',11,4)
par(mfrow=c(1,3))
for (i in (1:length(unique(geno)))){
	hist(acn.cs[[i]],main=unique(geno)[i])
	abline(v=c.score(t(acn[geno==unique(geno)[i],])),col='red',lty=2)
	obs=c.score(t(acn[geno==unique(geno)[i],]))
	p=length(acn.cs[[i]][acn.cs[[i]]<=obs])/length(acn.cs[[i]])
	legend('topright',legend=round(p,3))
	}

quartz('',11,5)
par(mfrow=c(1,3))
for (i in (1:length(acn.cup))){
	gplot(acn.cup[[i]],gmode='graph',mode='circle',vertex.col=apply(acn[geno==unique(geno)[i],],2,function(x) if (sum(x)!=0){1}else{'white'}),vertex.sides=100,vertex.border='grey',main=names(acn.cup)[i])
}
	


