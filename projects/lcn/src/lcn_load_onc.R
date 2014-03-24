###LCN: loading the onc data
###MKLau
###21Mar2014

source('~/projects/packages/cooc/src/cooc.R')
garden.data <- read.csv('~/projects/dissertation/projects/lcn/data/LCO_data_ONC_PIT.csv')
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
                                        #
if (all(table(onc[,1])==100)){print('Good to go!')}else{for (i in 1:1000){print('Warning: check input data!!!')}}
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
                                        #get genotype
onc.geno <- unlist(sapply(names(onc.q),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #Roughness in the Garden
rough <- read.csv('../../lcn/data/ONC_raw_roughness.csv')
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
load('../data/lcn_onc_ses.rda')
onc.ses <- os
if (all(names(onc.ses)==names(onc.q))){print('Good to go!')}else{print('Holy crap!')}
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
onc.rough <- avg.rough[match(ses.tree,r.tree)]
if (all(ses.tree==names(onc.rough))){print('Good to go!')}else{print('Holy Crap!')}
                                        #Microsat data from Nash
## gen.d <- read.csv(file='../../lcn/data/ONC_MSAT_datafromnash.csv')[,-1]
## gen.d[is.na(gen.d)] <- 0
## gen.d <- as.dist((gen.d))
                                        #RFLP distance values from Zink from Martinsen
rflp.d <- readLines('/Users/Aeolus/projects/dissertation/projects/acn/data/rflp_queller_goodnight.txt')
rflp.d <- lapply(rflp.d,strsplit,split='\t')
rflp.d <- lapply(rflp.d,function(x) x[[1]])
rflp.d[[61]] <- c(rflp.d[[61]],"")
rflp.d <- do.call(rbind,rflp.d)
rflp.n <- rflp.d[1,-1]
rflp.d <- rflp.d[-1,-1]
diag(rflp.d) <- 0
rflp.d <- matrix(as.numeric(rflp.d),nrow=nrow(rflp.d))
rownames(rflp.d) <- colnames(rflp.d) <- rflp.n
rflp.d <- as.dist(rflp.d)
                                        #community data
onc.com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
onc.R <- apply(sign(onc.com),1,sum)
onc.H <- diversity(onc.com)
onc.com.rel <- apply(onc.com,2,function(x) x/max(x))
onc.com.rel <- cbind(onc.com.rel,ds=rep(min(onc.com.rel[onc.com.rel!=0]),nrow(onc.com.rel)))
onc.com <- cbind(onc.com,ds=rep(min(onc.com[onc.com!=0]),nrow(onc.com)))
                                        #nmds and procrustes rotation
## onc.nms <- nmds.min(nmds(vegdist(onc.com.rel),2,2))
## onc.rot <- procrustes(dist(onc.rough),onc.nms,scale=TRUE)
## onc.rot <- t(onc.rot$rotation)
## onc.rot.sem <- onc.rot
## onc.rot <- onc.rot[,(1:2)[abs(cor(cbind(onc.rough,onc.rot))[1,2:3])==max(abs(cor(cbind(onc.rough,onc.rot))[1,2:3]))]]

### Renaming
oq <- onc.q
oc <- onc.com[,colnames(onc.com)!='ds']
os <- onc.ses
prb <- onc.rough
