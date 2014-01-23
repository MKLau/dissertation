#Reading Art's data files
#source('/Users/aeolus/Documents/R_Docs/EcoNet/EcoNet_03.R')
#source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/ind_net.R')
source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')

library(bipartite)
library(vegan)
library(sna)
library(gdata)

#Standard Error Function
se=function(x){sd(x)/sqrt(length(x))}

#Read Art's data files
read.art=function(file=""){
x<-read.xls(file,header=TRUE)
xspp<-x[,1]
x<-x[,-1]
x<-t(x)
colnames(x)<-xspp
return(x)
	}

my.split=function(x,lim=5){
y=unlist(strsplit(x,split=''))
x=character()
for (i in (1:min(c(lim,length(y))))){
	if (y[i]==' '|y[i]=='.'|y[i]=='X'){}else{
		x=paste(x,y[i],sep='')
		}
	}
as.character(x)
}

###NOTE!!! Fly (dance) in 2004 and 2005 was changed to Fly (Dolichopodidae) given the posibility of calling a Dolichopodidae, long legged fly, a "dance fly" and given that it was nearly un-detected in 2005 and not detected in 2004.

setwd("/Users/aeolus/Documents/R_Docs/Art_Keith/ArtONCArthropods")
ArtONC=list()
for (i in 1:3){
ArtONC[[i]]=read.art(dir()[i])
names(ArtONC)[i]=strsplit(dir()[i],split=" ")[[1]][1]
	}
names(ArtONC)
all(colnames(ArtONC[[1]])==colnames(ArtONC[[2]]))
all(colnames(ArtONC[[2]])==colnames(ArtONC[[3]]))

#bipartite graph perspective. These graphs show how each tree differs in its connections to the community of arthropods.
par(mfrow=c(1,3))
for (i in 1:length(ArtONC)){plotweb(ArtONC[[i]][order(apply(ArtONC[[i]],1,sum),decreasing=TRUE),order(apply(ArtONC[[i]],2,sum),decreasing=TRUE)],method='normal',text.rot=90)}

A=rbind(ArtONC[[1]],ArtONC[[2]],ArtONC[[3]])

#Cut out the species whose total abudances represent less than 1% of the total abundance across all species
barplot((apply(A,2,sum)/sum(apply(A,2,sum)))[order((apply(A,2,sum)/sum(apply(A,2,sum))),decreasing=TRUE)],las=2);abline(h=0.005)

A=rbind(ArtONC[[1]],ArtONC[[2]],ArtONC[[3]])
A=A[,(apply(A,2,sum)/sum(apply(A,2,sum)))>=0.005]

G=as.character(sapply(rownames(A),my.split))

A.list=list()
for (i in (1:length(unique(G)))){
	A.list[[i]]=A[G==unique(G)[i],]
	}
names(A.list)=unique(G)

A.net=list()

for (i in (1:length(A.list))){
A.net[[i]]=kendall.pairs(A.list[[i]],p.adj=FALSE,alpha=0.05,adj.method='fdr')
}

names(A.net)=names(A.list)

length(A.net)

quartz('ONC Arth. Networks',7,6.5)
par(mfrow=c(3,3))
for (i in (1:length(A.net))){
gplot(abs(A.net[[i]]),gmode='graph',mode='circle',vertex.col='black',edge.col='gray',vertex.sides=50,main=names(A.net)[i])	
}

##Network distance by hybrid index
#NOTE! There is no variation in hybrid index among these trees
hindex <- read.csv('/Users/Aeolus/Documents/garden_info/ONC_hybrid_index.csv')
hindex <- na.omit(hindex)
hi.geno <- as.character(gsub('-','',hindex[,1]))
hi <- as.numeric(hindex[,3])
geno.net <- as.character(names(A.net))
hi.net <- geno.net
geno.net[geno.net == 'Coal'] <- 'COAL1'
unique(geno.net)

for (i in (1:length(unique(geno.net)))){
hi.net[geno.net == unique(geno.net)[i]] <- hi[hi.geno == unique(geno.net)[i]]
}



##Figure for presentaiton
quartz('',22,22)
par(mfrow=c(3,3),mar=c(2.4,1.3,1.3,1.3),oma=c(0.1,0.1,0.1,0.1),bg='transparent',col.main='white',cex=2,mar=c(2,1,1,1))
net.list <- A.net
com.list <- A.list
                                        #names(net.list)=c('Excluded','Resistant','Susceptible') #rename the network graphs
net.list.reorder <- net.list
com.list.reorder <- com.list
for (i in (1:length(net.list))){
    v.col=apply(com.list.reorder[[i]],2,sum); v.col[v.col != 0] = 'lightblue' #color the present species
      v.col[v.col == 0] <- 'black' #color empty species
      gplot(abs(net.list.reorder[[i]]),gmode='graph',vertex.cex=3,vertex.sides=100,vertex.col=v.col,edge.lwd=0.35,edge.col=gray(0.9)[1],vertex.border='grey',mode='circle',displaylabels=FALSE,cex=2,main=names(net.list.reorder)[i]) #without titles
  }



#indicator species analyses
library(vegan)
library(labdsv)

G04=as.character(sapply(rownames(ArtONC[[1]]),my.split))
G05=as.character(sapply(rownames(ArtONC[[2]]),my.split))
G06=as.character(sapply(rownames(ArtONC[[3]]),my.split))

summary(indval(ArtONC[[1]],G04,numitr=5000))
G04[order(G04)]
summary(indval(ArtONC[[2]],G05,numitr=5000))
G04[order(G05)]
summary(indval(ArtONC[[3]],G06,numitr=5000))
G04[order(G06)]

par(mfrow=c(1,3))
plot(ArtONC[[1]][,colnames(ArtONC[[1]])=='Leaf miner']~factor(G04),las=2,ylab='Leaf miner',xlab='')
plot(ArtONC[[1]][,colnames(ArtONC[[1]])=='Stem Borer']~factor(G04),las=2,ylab='Stem Borer',xlab='')
plot(ArtONC[[1]][,colnames(ArtONC[[1]])=='T. Dip (tiny)']~factor(G04),las=2,ylab='T. Dip (tiny)',xlab='')

par(mfrow=c(2,4))
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Leafhopper (brwn mot)']~factor(G05),las=2,ylab='Leafhopper (brwn mot)',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Beetle (unknown)']~factor(G05),las=2,ylab='Beetle (unknown)',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Coenagrionidae']~factor(G05),las=2,ylab='Coenagrionidae',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Leafhopper (fsh eye)']~factor(G05),las=2,ylab='Leafhopper (fsh eye)',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Parasitoid']~factor(G05),las=2,ylab='Parasitoid',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Calophorid']~factor(G05),las=2,ylab='Calophorid',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Leafhopper (gry str)']~factor(G05),las=2,ylab='Leafhopper (gry str)',xlab='')
plot(ArtONC[[2]][,colnames(ArtONC[[2]])=='Thrip (black)']~factor(G05),las=2,ylab='Thrip (black)',xlab='')

par(mfrow=c(1,2))
plot(ArtONC[[3]][,colnames(ArtONC[[3]])=='Beetle (unknown)']~factor(G06),las=2,xlab='',ylab='Beetle (unknown)')
plot(ArtONC[[3]][,colnames(ArtONC[[3]])=='Salticidae (gry/blk)']~factor(G06),las=2,xlab='',ylab='Salticidae (gry/blk)')


#Using across years as replicates
#NOTE: this didn't work, most likely because of an issue of low power.
A=rbind(ArtONC[[1]],ArtONC[[2]],ArtONC[[3]])
A=A[,(apply(A,2,sum)/sum(apply(A,2,sum)))>=0.005]

A.tree=list()
for (i in (1:length(unique(rownames(A))))){
	A.tree[[i]]=A[rownames(A)==unique(rownames(A))[i],]
}

names(A.tree)=unique(rownames(A))
G.tree=as.character(sapply(names(A.tree),my.split))

Atree.nets=list()
for (i in (1:length(A.tree))){
	Atree.nets[[i]]=kendall.pairs(A.tree[[i]],p.adj=FALSE,alpha=0.1)
	Atree.nets[[i]][is.na(Atree.nets[[i]])]=0
}

for (i in (1:length(Atree.nets))){
	gplot(abs(Atree.nets[[i]]))
	locator(1)
}
	


