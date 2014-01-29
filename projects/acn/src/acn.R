###Final set of analyses for the ACN paper
###29 Jan 2014
###To be run up on Hoth

library(vegan)
library(bipartite)
library(pbapply)

###Sourcing ComGenR on hoth
oldwd <- getwd()
setwd('../../../../ComGenR/R/')
cgn.list <- (sapply(dir(),grepl,pattern='~')|sapply(dir(),grepl,pattern='\\#'))==FALSE
sapply(dir()[cgn.list],source)
setwd(oldwd)

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
liv.geno <- as.character(unlist(lapply(split(geno$liv,split(tree,pit[,1])$liv),function(x) x[1])))
sen.geno <- as.character(unlist(lapply(split(geno$sen,split(tree,pit[,1])$sen),function(x) x[1])))
liv.com <- do.call(rbind,lapply(split(liv,split(tree,pit[,1])$liv),function(x) apply(x,2,sum)))
sen.com <- do.call(rbind,lapply(split(sen,split(tree,pit[,1])$sen),function(x) apply(x,2,sum)))
tree <- pit$tree
tree.geno <- as.character(unlist(lapply(split(unlist(geno),tree),function(x) x[1])))
leaf.type <- as.character(pit[,1])

## Genetic effect on P. betae
                                        #total PB across live and senescent leaves
pb.total <- tapply(pit$pb,tree,sum)
summary(aov(pb.total~tree.geno))

## Genotype effect on composition
adonis(liv.com~liv.geno)
liv.com.rel <- apply(liv.com,2,function(x) x/max(x))
liv.com.rel[is.na(liv.com.rel)] <- 0
adonis(liv.com~liv.geno)
adonis(liv.com.rel~liv.geno)
sen.com. <- cbind(sen.com,ds=rep(1,nrow(sen.com)))
sen.com.rel <- apply(sen.com.,2,function(x) x/max(x))
sen.com.rel[,ncol(sen.com.rel)] <- min(sen.com.rel[sen.com.rel!=0])/10
adonis(sen.com.~sen.geno)
adonis(sen.com.rel~sen.geno)

## Genetic effect on SES
liv.ses <- cnm.test(liv.com,nits=1000)
sen.ses <- cnm.test(sen.com,nits=1000)

## Co-occurrence network structure


## Nestedness for live and senesced
liv.gnet <- mean.g(liv.com,liv.geno)
sen.gnet <- mean.g(sen.com,sen.geno)
par(mfrow=c(1,2))
cgPlotweb(liv.com,liv.geno)
cgPlotweb(sen.com,sen.geno)
oecosimu(liv.gnet,method='r0',nestfun=nestedtemp,nits=1000,burn=50)
oecosimu(sen.gnet,method='r0',nestfun=nestedtemp,nits=1000,burn=50)
