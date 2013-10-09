###Wild stand co-occurrence analyses

###Tree data
x <- read.csv('~/projects/dissertation/projects/lichen_coo/data/lco_Apr2012.csv')
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
                                        #condense species
lec.spp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x[,-ncol(x)],lec=lec.spp,phy=phy.spp,NOTES=x[,ncol(x)])
x <- x[,-c(6,8,10,18,13,14,15,16)]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
x.q <- split(x,quads)
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
stand <- x[,c(-1:-4,-ncol(x))]
stand.null <- nullCom(stand)
stand.obs <- cscore(stand)
stand.cs <- unlist(lapply(stand.null,cscore))
stand.ses <- (mean(stand.cs)-stand.obs)/sd(stand.cs)
stand.p <- length(stand.cs[stand.cs<stand.obs])/length(stand.obs)
stand.out <- cbind(stand.ses=stand.ses,stand.p=stand.p)
write.csv(stand.out,file='~/projects/dissertation/projects/lichen_coo/results/wild_coo_results.csv')
