###ACN testing for nestedness
###30 Jan 2014

library(vegan)
library(bipartite)
library(pbapply)

###Sourcing ComGenR on hoth
oldwd <- getwd()
setwd('../../../../ComGenR/R/')
cgn.list <- (sapply(dir(),grepl,pattern='~')|sapply(dir(),grepl,pattern='\\#'))==FALSE
sapply(dir()[cgn.list],source)
setwd(oldwd)

###Data
## Load data
pit <- read.csv('~/projects/dissertation/projects/acn/data/arth_cooc_PIT_Lau.csv')
                                        #genotype remove 1007
pit <- pit[pit$geno!='1007',]
pit$geno <- factor(as.character(pit$geno))
pit[is.na(pit)] <- 0
                                        #remove trailing 0
pit$tree <- sub('\\.0','\\.',pit$tree)
                                        #data verification of trees and genotypes
pg <- read.csv('~/projects/dissertation/docs/garden_information/PIT_garden_tree_information.csv')
pg <- pg[pg$Row!='DEAD',]
pg.tree <- as.character(pg$Row)
pg.tree <- sub('-','.',pg.tree)
pg.tree <- sub('NP','np',pg.tree)

                                        #pemphigus mergers
pit$pb.upper <- pit$pb.upper + pit$pb.woody
pb <- pit$pb.upper + pit$pb.lower
pit <- cbind(pit,pb=pb)
pit$pb.pred <- pit$pb.pred + pit$pb.hole + pit$pb.woody.pred
pit <- pit[,colnames(pit)!='pb.woody'&colnames(pit)!='pb.woody.pred'&colnames(pit)!='pb.hole'&colnames(pit)!='mite'&colnames(pit)!='pb.upper'&colnames(pit)!='pb.lower']
                                        #separate live and senescing leaves
pit.com <- pit[,-1:-7] #this removes tree info as well as fungus
liv <- pit.com[pit[,1]=='live',]
sen <- pit.com[pit[,1]=='sen',]
geno <- split(pit[,3],pit[,1])
tree <- pit$tree
liv.geno <- as.character(unlist(lapply(split(geno$liv,split(tree,pit[,1])$liv),function(x) x[1])))
sen.geno <- as.character(unlist(lapply(split(geno$sen,split(tree,pit[,1])$sen),function(x) x[1])))
liv.com <- do.call(rbind,lapply(split(liv,split(tree,pit[,1])$liv),function(x) apply(x,2,sum)))
sen.com <- do.call(rbind,lapply(split(sen,split(tree,pit[,1])$sen),function(x) apply(x,2,sum)))
tree.geno <- as.character(unlist(lapply(split(unlist(geno),tree),function(x) x[1])))
leaf.type <- as.character(pit[,1])
tree.leaf <- unlist(lapply(split(leaf.type,paste(leaf.type,tree)),function(x) x[1]))
n.leaves <- unlist(lapply(split(pit.com,paste(leaf.type,tree)),function(x) nrow(x)))
liv.pc <- liv.com
for (i in 1:nrow(liv.com)){liv.pc[i,] <- liv.com[i,]/split(n.leaves,tree.leaf)[[1]][i]}
sen.pc <- sen.com
for (i in 1:nrow(sen.com)){sen.pc[i,] <- sen.com[i,]/split(n.leaves,tree.leaf)[[2]][i]}
pit.trees <- split(pit.com,paste(leaf.type,tree))

###Nestedness
liv.gnet <- mean.g(liv.com,liv.geno)
sen.gnet <- mean.g(sen.com,sen.geno)
                                        #live
liv.nest <- list()
liv.nest[[1]] <- oecosimu(liv.gnet,method='r00',nestfun=nestedtemp,nits=5000,burn=50)
liv.nest[[2]] <- oecosimu(liv.gnet,method='r0',nestfun=nestedtemp,nits=5000,burn=50)
liv.nest[[3]] <- oecosimu(liv.gnet,method='c0',nestfun=nestedtemp,nits=5000,burn=50)
liv.nest[[4]] <- oecosimu(liv.gnet,method='r1',nestfun=nestedtemp,nits=5000,burn=50)
names(liv.nest) <- c('r00','r0','c0','r1')
liv.nest <- do.call(rbind,lapply(liv.nest,function(x) as.numeric(unlist(x$oecosimu)[c(6,1,3)])))
colnames(liv.nest) <- c('temp','z','P')
                                        #write to disk
write.csv(liv.nest,file='../results/nest_liv.csv')
                                        #sen
sen.nest <- list()
sen.nest[[1]] <- oecosimu(sen.gnet,method='r00',nestfun=nestedtemp,nits=5000,burn=50)
sen.nest[[2]] <- oecosimu(sen.gnet,method='r0',nestfun=nestedtemp,nits=5000,burn=50)
sen.nest[[3]] <- oecosimu(sen.gnet,method='c0',nestfun=nestedtemp,nits=5000,burn=50)
sen.nest[[4]] <- oecosimu(sen.gnet,method='r1',nestfun=nestedtemp,nits=5000,burn=50)
names(sen.nest) <- c('r00','r0','c0','r1')
sen.nest <- do.call(rbind,lapply(sen.nest,function(x) as.numeric(unlist(x$oecosimu)[c(6,1,3)])))
colnames(sen.nest) <- c('temp','z','P')
                                        #write to disk
write.csv(sen.nest,file='../results/nest_sen.csv')
