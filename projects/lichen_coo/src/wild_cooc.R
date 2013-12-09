###Wild stand co-occurrence analyses

###Tree data
x <- read.csv('~/projects/dissertation/projects/lichen_coo/data/lco_Apr2012.csv')
                                        #remove notes
x <- x[,colnames(x)!='NOTES.']
                                        #
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
                                        #condense species
lec.spp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x,lec=lec.spp,phy=phy.spp)
x <- x[,-c(6,8,10,18,13,14,15,16)]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
x.q <- split(x,quads)
x.t <- split(x,as.character(x$tree))
                                        #environmental data from lamit
env <- read.csv('~/projects/dissertation/projects/lichen_coo/data/Uinta2012_all_data_from_Lamit.csv')
env <- env[is.na(env$Pct.Roughness)==FALSE,]
env[,1] <- sub('\\?','',sub('\\.0','\\.',sub('\\_','\\.',sub('\\-','\\.',tolower(as.character(env[,1]))))))
env[env[,1]=='ll.6_(znu.29)',1] <- 'll.6'
env[env[,1]=='gnu.85.1ftaway',1] <- 'gnu.85'
env$Quad.Loc <- as.character(sapply(as.character(env$Quad.Loc),function(x) unlist(strsplit(x,split='_'))[2]))
env$Quad.Loc <- sub('\\-','\\.',env$Quad.Loc)
env$Quad.Loc <- paste('n',env$Quad.Loc,sep='')
                                        #remove southern aspect
env <- env[env$Aspect!='South',]
env.tid <- paste(env$Tree.ID,env$Quad.Loc)
                                        #check that the datasets are compatible
all(names(x.q)%in%env.tid) 
                                        #match observations
all(env.tid[match(names(x.q),env.tid)] == names(x.q))
                                        #delimit to co-occurrence dataset and match
env <- env[match(names(x.q),env.tid),]
                                        #
library(vegan)
###Co-occurrence C-score
##Whole Stand
source('~/projects/dissertation/projects/lichen_coo/src/seenetR.R')
                                        #
## stand <- x[,c(-1:-4,-ncol(x))]
## stand.null <- nullCom(stand)
## stand.obs <- cscore(stand)
## stand.cs <- unlist(lapply(stand.null,cscore))
## stand.ses <- (stand.obs-mean(stand.cs))/sd(stand.cs)
## stand.p <- length(stand.cs[stand.cs<stand.obs])/length(stand.obs)
## stand.out <- cbind(stand.ses=stand.ses,stand.p=stand.p)
## write.csv(stand.out,file='~/projects/dissertation/projects/lichen_coo/results/wild_coo_results.csv')
#Tree scale
nits <- 5000
x.t. <- lapply(x.t,function(x) x[,-1:-4])
tree.ocs <- unlist(lapply(x.t.,cscore))
tree.ncs <- list() #null communities
tree.ses <- tree.pval <- numeric(length(x.t.))
for (i in 1:length(x.t)){
  if (sum(x.t.[[i]])<1){tree.nc[[i]] <- NA}else{
    tree.ncs[[i]] <- unlist(lapply(nullCom(x.t.[[i]],nits=nits),cscore))
  }
  tree.ses[i] <- (tree.ocs[i] - mean(tree.ncs[[i]]))/sd(tree.ncs[[i]])
  tree.pval[i] <- length(tree.ncs[[i]][tree.ncs[[i]]<=tree.ocs[i]])/length(tree.ncs[[i]])
  print(i)
}
                                        #output
dput(tree.ocs,'../data/wild_tree_ocs.Rdata')
dput(tree.ncs,'../data/wild_tree_ncs.Rdata')
dput(tree.ses,'../data/wild_tree_ses.Rdata')
dput(tree.pval,'../data/wild_tree_pval.Rdata')
