###Arthropod Co-occurrence Patterns
###Gina and Art's datasets

library(xlsx)
library(vegan)
library(ecodist)
library(rgl)
library(mvabund)
library(sna)
#source(ComGenR.R

###Gina's data
                                        #using 2003 data only
gina <- read.xlsx('~/data/2000-2003 Garden Data from Gina.xls',4)
                                        #remove backcrosses
xtype <- substr(gina[,1],1,2)
gina <- gina[,-1:-2]
                                        #make NAs zeros
gina[is.na(gina)] <- 0
                                        #
###Composition
com <- gina[,apply(gina,2,sum)!=0]
adonis(com~xtype)

###Unipartite network
net <- CoNetwork(com[,apply(com,2,sum)!=0])
n.com <- nullCom(com[,apply(com,2,sum)!=0])
n.c <- unlist(lapply(n.com,C.score))
o.c <- C.score(com[,apply(com,2,sum)!=0])
p.val <- length(n.c[n.c<=o.c])/length(n.c)
ses <- (o.c-mean(n.c))/sd(n.c)
                                        #
v.cex <- apply(com[,apply(com,2,sum)!=0],2,sum) #scaling node size by the log of species frequencies
v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*3)+0.1
e.col <- net
e.col[net>0.5] <- 'red'
e.col[net<0.5] <- 'black'
e.col[net==0.5] <- 'grey'
gplot(abs(net),displaylabels=FALSE,gmode='graph',pad=1.5,
      edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',
      edge.col=e.col)

###Bipartite Network
xtype. <- xtype[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE)]
bpn <- com[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE),order(apply(com,2,function(x) sum(sign(x))),decreasing=TRUE)]
rownames(bpn) <- paste(xtype.,(1:nrow(bpn)),sep='_')
plotweb(bpn,method='normal',text.rot=90,col.low=as.numeric(factor(xtype.)))
nest.test <- oecosimu(bpn, nestfun="nestedtemp", "r1",nsimul = 999,burnin=50)
nest.test <- oecosimu(bpn, nestfun="nestedtemp", "r00",nsimul = 999,burnin=50)

###Modularity
bpn.mods <- computeModules(bpn[,apply(bpn,2,function(x) sum(sign(x)))>1])
plotModuleWeb(bpn.mods)

###Plot
par(mfcol=c(2,1))
gplot(abs(net),displaylabels=FALSE,gmode='graph',pad=1.5,
      edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',
      edge.col=e.col)
plotweb(bpn,method='normal',text.rot=90,col.low=as.numeric(factor(xtype.)))

###Art's Data
art <- list()
art[[1]] <- read.csv('~/data/onc_arth_04.csv')[,-1]
art[[2]] <- read.csv('~/data/onc_arth_05.csv')[,-1]
art[[3]] <- read.csv('~/data/onc_arth_06.csv')[,-1]
geno <- read.csv('~/data/onc_arth_04.csv')[,1]
year <- c(rep('2004',nrow(art[[1]])),rep('2005',nrow(art[[2]])),rep('2006',nrow(art[[3]])))
sp.rm <- apply(do.call(rbind,lapply(art,function(x) apply(x,2,function(x) sum(sign(x))))),2,function(x) all(x!=0))
art <- do.call(rbind,art)
art <- art[,sp.rm]
art <- split(art,year)
                                        #remove species that are not present in all three years
                                        #analysis
library(igraph)
par(mfrow=c(1,3))
art.net <- lapply(art,CoNetwork)
                                        #test genotype effects
art.rel <- lapply(art,function(x) apply(x,2,function(x) x/max(x)))
art.perm <- lapply(art,function(x,g) adonis(x~g)$aov.tab,g=geno)
artrel.perm <- lapply(art.rel,function(x,g) adonis(x~g)$aov.tab,g=geno)
do.call(rbind,art.perm)
do.call(rbind,artrel.perm)
