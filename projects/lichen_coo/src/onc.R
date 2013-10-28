###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????
rm(list=ls())
library(sna)
source('./helper_funcs.R')

###Garden Analysis
garden.data <- read.csv('~/projects/dissertation/projects/lichen_coo/data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'

                                        #separate onc
onc <- garden.data[garden=='onc',]
					#tree overlap between years
unique(onc$Tree[onc$Year=='2010']) %in% unique(onc$Tree[onc$Year=='2011'])
unique(onc$Tree[onc$Year=='2011']) %in% unique(onc$Tree[onc$Year=='2010'])
                                        #onc <- onc[onc$Year=='2011',]
###Merge species groups
Phy <- apply(onc[,12:14],1,sum) #make phy out of pmel,pads,pund
onc. <- onc[,-12:-14]
onc <- data.frame(onc.,Phy)
###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env);adonis(com~height);detach(env)
library(pbapply)
source('./seenetR.R')
###Co-occurrences
##stand level
                                        #onc network
onc.cn <- co.net(onc[,7:ncol(onc)])
onc.dn <- dep.net(onc[,7:ncol(onc)])
onc.graph <- onc.dn[apply(onc.dn,1,sum)!=0,apply(onc.dn,2,sum)!=0]
par(mfrow=c(1,2))
my.gplot(onc.dn,v.cex=((apply(com,2,sum)/max(apply(com,2,sum)))+0.5))
onc.deg <- degree(onc.dn)
names(onc.deg) <- rownames(onc.cn)
barplot(onc.deg,ylab='Centrality')
                                        #co-occurrence
stand.null <- unlist(dget(file='../data/onc_stand_null.Rdata'))
stand.ses <- (cscore(onc[,7:ncol(onc)]) - mean(stand.null)) / sd(stand.null)
stand.ses.p <- length(stand.null[stand.null<=stand.ses])/length(stand.null)
c(stand.ses,stand.ses.p)
##tree level 
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,3],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
obs.cs <- unlist(lapply(onc.q,cscore))
onc.ses <- dget('../data/onc_tree_ses.Rdata')
onc.p <- dget('../data/onc_tree_pval.Rdata')
ses.zero.p <- FALSE
ses.zero.sd2 <- TRUE
if (ses.zero.p){onc.ses[onc.p>0.05] <- 0}else{}
if (ses.zero.sd2){onc.ses[abs(onc.ses) < 2] <- 0}else{}
onc.tn <- lapply(onc.q,dep.net) #tree level networks
names(onc.ses) <- names(onc.q)
onc.ses[is.na(onc.ses)] <- 0
onc.centrality <- unlist(lapply(onc.tn,function(x) centralization(x,FUN='degree')))
onc.deg <- unlist(lapply(onc.tn,function(x) length(x[x!=0])))
###Roughness in the Garden
rough <- read.csv('../data/ONC_raw_roughness.csv')
rough <- rough[as.character(rough[,1])!="",1:5]
                                        #isolate north quadrats
rough <- rough[sapply(rough[,3],function(x) substr(x,1,1)=='N'),]
                                        #average roughness
avg.rough <- tapply(rough[,5],rough[,1],mean)
r.tree <- names(avg.rough)
r.tree <- sub('-','\\.',r.tree)
r.tree <- sub('\\.0','\\.',r.tree)
names(avg.rough) <- r.tree
                                        #match roughness to to ses values
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
avg.rough <- avg.rough[match(ses.tree,r.tree)]
all(ses.tree==names(avg.rough))

###Genotype
genotype <- as.character(sapply(names(onc.q),function(x) unlist(strsplit(x,split=' '))[3]))
                                        #remove outlier trees
rm.n1.31 <- (names(onc.ses)!='N1.31 2011 1008')
rm.rl6 <- (genotype!='RL6')
                                        #roughness by genotype
summary(aov(avg.rough~genotype))
summary(aov(avg.rough[rm.n1.31&rm.rl6]~genotype[rm.n1.31&rm.rl6]))
                                        #centrality
plot(sqrt(onc.centrality)~factor(genotype))
summary(aov(onc.ses[names(onc.ses)!='N1.31 2011 1008']~genotype[names(onc.ses)!='N1.31 2011 1008']))
plot(onc.ses[names(onc.ses)!='N1.31 2011 1008']~factor(genotype[names(onc.ses)!='N1.31 2011 1008']))
                                        #
geno.rm <- c('H10','RL6','T6','WC5')
rm.n1.31 <- (names(onc.ses)!='N1.31 2011 1008')
rm.rl6 <- (genotype!='RL6')
outlier.rm <- ((names(onc.ses)!='N1.31 2011 1008')&(genotype%in%geno.rm==FALSE))
summary(aov(onc.ses[outlier.rm]~genotype[outlier.rm]))
plot(onc.ses[outlier.rm]~factor(genotype[outlier.rm]))
                                        #
summary(aov(sqrt(onc.centrality[outlier.rm])~genotype[outlier.rm]))
plot(onc.centrality[outlier.rm]~factor(genotype[outlier.rm]))
                                        #
summary(aov(log(onc.deg+0.001)~genotype))
plot(onc.deg~factor(genotype))
summary(aov(log(onc.deg[outlier.rm]+0.001)~genotype[outlier.rm]))
plot(onc.deg[outlier.rm]~factor(genotype[outlier.rm]))
###Removed N1.31 and RL6
                                        #roughness ~ genotype 
summary(aov(avg.rough[rm.n1.31&rm.rl6]~genotype[rm.n1.31&rm.rl6]))
plot(avg.rough[rm.n1.31&rm.rl6]~factor(genotype[rm.n1.31&rm.rl6]))
                                        #ses and centrality correlation
summary(lm(onc.ses[rm.n1.31&rm.rl6]~onc.deg[rm.n1.31&rm.rl6]))
plot(onc.ses[rm.n1.31&rm.rl6]~onc.deg[rm.n1.31&rm.rl6])
                                        #
summary(aov((log(abs(onc.ses[rm.n1.31&rm.rl6])+0.1))~genotype[rm.n1.31&rm.rl6]))
plot(onc.ses[rm.n1.31&rm.rl6]~factor(genotype[rm.n1.31&rm.rl6]))
summary(aov((log(abs(onc.ses[rm.n1.31&rm.rl6])+0.1))~avg.rough[rm.n1.31&rm.rl6]))
                                        #
summary(aov(sqrt(onc.centrality[rm.n1.31&rm.rl6])~genotype[rm.n1.31&rm.rl6]))
plot(onc.centrality[rm.n1.31&rm.rl6]~factor(genotype[rm.n1.31&rm.rl6]))
summary(aov(sqrt(onc.centrality[rm.n1.31&rm.rl6])~avg.rough[rm.n1.31&rm.rl6]))
plot(onc.centrality[rm.n1.31&rm.rl6]~factor(avg.rough[rm.n1.31&rm.rl6]))
summary(aov((onc.centrality[rm.n1.31&rm.rl6])~avg.rough[rm.n1.31&rm.rl6]))
                                        #
summary(aov(log(onc.deg[rm.n1.31&rm.rl6]+0.001)~genotype[rm.n1.31&rm.rl6]))
plot(onc.deg[rm.n1.31&rm.rl6]~factor(genotype[rm.n1.31&rm.rl6]))
summary(lm((onc.deg[rm.n1.31&rm.rl6])^2~avg.rough[rm.n1.31&rm.rl6]))
plot(onc.deg[rm.n1.31&rm.rl6]~avg.rough[rm.n1.31&rm.rl6])
                                        #deg~rough
summary(aov(sqrt(onc.deg[rm.n1.31&rm.rl6])~avg.rough[rm.n1.31&rm.rl6]))                                        
                                        #network similarity
net.dist <- matrix(0,nrow=length(onc.tn[rm.n1.31&rm.rl6]),ncol=length(onc.tn[rm.n1.31&rm.rl6]))
for (i in 1:nrow(net.dist)){
  for (j in 1:ncol(net.dist)){
    x <- onc.tn[rm.n1.31&rm.rl6][[i]]
    y <- onc.tn[rm.n1.31&rm.rl6][[j]]
    net.dist[i,j] <- sqrt(sum((x - y)^2))
  }
}
net.d <- as.dist((net.dist))
adonis(net.d~factor(genotype[rm.n1.31&rm.rl6]),perm=5000)
adonis(net.d~avg.rough[rm.n1.31&rm.rl6],perm=5000)
                                        #composition
adonis(com.rel[rm.n1.31&rm.rl6,]~avg.rough[rm.n1.31&rm.rl6],perm=5000)
                                        #look at compositional effect of genotype
com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
com. <- cbind(com,ds=rep(1,nrow(com)))
com.rel <- apply(com.,2,function(x) x/max(x))
com.rel[,ncol(com.rel)] <- com.rel[,ncol(com.rel)] * 0.001
adonis(com.~factor(genotype))
adonis(com.~avg.rough)
adonis(com.rel~factor(genotype))
adonis(com.rel~avg.rough)

###SEM
###Can we get genotype averages?

###NEED TO RECONCILE TAXONOMIC DIFFERENCES
