###Arthropod Co-occurrence Patterns
###Gina and Art's datasets

library(ComGenR)

###Gina's data
                                        #using 2003 data only
gina <- read.csv('~/data/2000-2003 Garden Data from Gina.csv')
xtype <- substr(gina[,1],1,2)
gina <- gina[,-1:-2]
gina[is.na(gina)] <- 0 
                                        #
###Composition
com <- gina[,apply(gina,2,sum)!=0]
adonis(com~xtype)

###Unipartite network
cnm.results <- cnm.test(com,nits=1000)
net <- CoNetwork(com[,apply(com,2,sum)!=0],threshold=0)
mgp(net,com)

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

###Art's Data
###Control trees from Art's PB experiment can be used for non-Art-first-author work.
###Analyze both for Art's pub.

setwd('~/projects/pb_removal_nets/')
art <- list(read.csv('./data/keith_pb_removal_2008.csv'),read.csv('./data/keith_pb_removal_2009.csv'))
art[[1]][is.na(art[[1]])] <- 0
art[[2]][is.na(art[[2]])] <- 0
art <- lapply(art,function(x) split(x,(1:nrow(x))%%2))
art <- list(art[[1]][[1]],art[[1]][[2]],art[[2]][[1]],art[[2]][[2]])
names(art) <- c('2008c','2008x','2009c','2009x')
                                        #permanova
art.perm <- lapply(art,function(x) adonis(x[,-1:-2]~x[,1])$aov.tab)
do.call(rbind,art.perm)
artrel.perm <- lapply(art,function(x) adonis(apply(x[,-1:-2],2,function(x) x/max(x))[,apply(x[,-1:-2],2,sum)!=0]~x[,1])$aov.tab)
do.call(rbind,artrel.perm)

###
###For the EEN paper use 2008 control, as it doesn't have the issue of zero sum trees
###
                                        #unipartite networks
art.nets <- lapply(lapply(art,function(x) x[,-1:-2]),CoNetwork,threshold=3)
unlist(lapply(art.nets,function(x) length(x[x!=0])))
art.minnet <- min.net(art.nets[[1]],art[[1]][,-1:-2])
mgp(art.nets[[1]],art[[1]][,-1:-2],displaylabels=FALSE)
mgp(art.minnet[[1]],art.minnet[[2]],displaylabels=TRUE)


mark

##Test for composition


##Mantel of species distances between two networks




## art <- list()
## art[[1]] <- read.csv('~/data/onc_arth_04.csv')[,-1]
## art[[2]] <- read.csv('~/data/onc_arth_05.csv')[,-1]
## art[[3]] <- read.csv('~/data/onc_arth_06.csv')[,-1]
## geno <- read.csv('~/data/onc_arth_04.csv')[,1]
## year <- c(rep('2004',nrow(art[[1]])),rep('2005',nrow(art[[2]])),rep('2006',nrow(art[[3]])))
## sp.rm <- apply(do.call(rbind,lapply(art,function(x) apply(x,2,function(x) sum(sign(x))))),2,function(x) all(x!=0))
## art <- do.call(rbind,art)
## art <- art[,sp.rm]
## art <- split(art,year)
##                                         #remove species that are not present in all three years
##                                         #analysis
## library(igraph)
## par(mfrow=c(1,3))
## art.net <- lapply(art,CoNetwork)
##                                         #test genotype effects
## art.rel <- lapply(art,function(x) apply(x,2,function(x) x/max(x)))
## art.perm <- lapply(art,function(x,g) adonis(x~g)$aov.tab,g=geno)
## artrel.perm <- lapply(art.rel,function(x,g) adonis(x~g)$aov.tab,g=geno)
## do.call(rbind,art.perm)
## do.call(rbind,artrel.perm)
