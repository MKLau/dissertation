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
                                        #remove genotype RL6 and N1.31
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Tree!='N1.31',]
                                        #separate onc
garden.data[,1] <- as.character(garden.data[,1])
g1 <- substr(garden.data[,1],2,2)
g1[g1!='P'] <- 'onc'
onc <- garden.data[g1=='onc',]
					#tree overlap between years
unique(onc$Tree[onc$Year=='2010']) %in% unique(onc$Tree[onc$Year=='2011'])
unique(onc$Tree[onc$Year=='2011']) %in% unique(onc$Tree[onc$Year=='2010'])
                                        #onc <- onc[onc$Year=='2011',]
###Merge species groups
Phy <- apply(onc[,12:14],1,sum) #make phy out of pmel,pads,pund
onc. <- onc[,-12:-14]
onc <- data.frame(onc.,Phy)
if (all(table(onc[,1])==100)){}else{for (i in 1:1000){print('Warning: check input data!!!')}}
###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env);adonis(com~height);detach(env)

###modeling
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
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
obs.cs <- unlist(lapply(onc.q,cscore))
                                        #load ses values
onc.ses <- dget('../data/onc_tree_ses.Rdata')
onc.p <- dget('../data/onc_tree_pval.Rdata')
ses.zero.p <- TRUE;ses.zero.sd2 <- FALSE
if (ses.zero.p){onc.ses[onc.p>0.05] <- 0}else{}
if (ses.zero.sd2){onc.ses[abs(onc.ses) < 2] <- 0}else{}
onc.tn <- lapply(onc.q,dep.net) #tree level networks
names(onc.ses) <- names(onc.q)[names(onc.q)!='']
onc.ses[is.na(onc.ses)] <- 0
onc.ses <- onc.ses[is.na(names(onc.ses))!=TRUE]
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
genotype <- as.character(sapply(names(onc.q),function(x) unlist(strsplit(x,split=' '))[2]))
###Analyses
onc.ses
