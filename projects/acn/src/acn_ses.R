###ACN ses values for trees
###30Jan2014

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

###SES values
pit.trees <- split(pit.com,paste(leaf.type,tree))
ses.trees <- lapply(pit.trees,cnm.test)
ses.trees <- do.call(rbind,ses.trees)
write.csv(ses.trees,'../data/acn_ses.csv')

###Nestedness values
liv.nest <- read.csv('../results/nest_liv.csv')
sen.nest <- read.csv('../results/nest_sen.csv')
