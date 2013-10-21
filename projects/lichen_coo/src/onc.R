###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????
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

###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env);adonis(com~height);detach(env)

###SES values
library(pbapply)
source('./seenetR.R')

##stand level
                                        #onc network
onc.cn <- co.net(onc[,7:15])
onc.dn <- dep.net(onc[,7:15])
onc.graph <- onc.dn[apply(onc.dn,1,sum)!=0,apply(onc.dn,2,sum)!=0]
my.gplot(onc.dn)
onc.deg <- degree(onc.dn)
names(onc.deg) <- rownames(onc.cn)
barplot(onc.deg,ylab='Centrality')
                                        #co-occurrence
stand.null <- unlist(dget(file='../data/onc_stand_null.Rdata'))
stand.ses <- (cscore(onc[,7:15]) - mean(stand.null)) / sd(stand.null)
stand.ses.p <- length(stand.null[stand.null<=stand.ses])/length(stand.null)
c(stand.ses,stand.ses.p)
##tree level 
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,3],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:15])
obs.cs <- unlist(lapply(onc.q,cscore))
onc.ses <- dget('../data/onc_tree_ses_2011.Rdata')
names(onc.ses) <- names(onc.q)
     # onc.sim <- pblapply(onc.q,function(x) if (sum(sign(apply(x,2,sum)))>1){nullCom(x)}else{NA})
     # onc.cs <- pblapply(onc.sim,function(x) if (any(is.na(x[[1]]))){NA}else{lapply(x,cscore)})
     # onc.cs <- pblapply(onc.cs,unlist)
     # onc.ses <- obs.cs*0
     # for (i in 1:length(onc.ses)){
     #   onc.ses[i] <- (obs.cs[i] - mean(onc.cs[[i]])) / sd(onc.cs[[i]])
     # }

###Genotype
genotype <- as.character(sapply(names(onc.q),function(x) unlist(strsplit(x,split=' '))[3]))
onc.ses[is.na(onc.ses)] <- 0
summary(aov(onc.ses~genotype))
plot(onc.ses~factor(genotype))

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
                                        #match to ses values
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
avg.rough <- avg.rough[match(ses.tree,r.tree)]
all(ses.tree==names(avg.rough))
                                        #
summary(lm(onc.ses~avg.rough))

                                        #look at compositional effect of genotype
com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
com. <- cbind(com,ds=rep(1,nrow(com)))
adonis(com.~factor(genotype))
adonis(com.~avg.rough)

###NEED TO RECONCILE TAXONOMIC DIFFERENCES
