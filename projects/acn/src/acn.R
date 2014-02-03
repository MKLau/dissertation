###Final set of analyses for the ACN paper
###29 Jan 2014
###To be run up on Hoth

###Re-do labeling so that any thing at the:
##leaf scale has no suffix
##tree (add up the leaves for each tree) scale is labeled t
##int. (values summed up for each tree) is labled i
rm(list=ls())
library(vegan)
library(bipartite)
library(pbapply)
library(ComGenR)
                                        #reml
                                        #using notes from the following link
                                        #http://www.stat.wisc.edu/~ane/st572/notes/lec21.pdf
source('~/projects/packages/ComGenR_development/src/cgREML.R')
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
geno <- pit[,3]
tree <- pit[,2]
type <- pit[,1]
obs <- split(pit.com,paste(tree,type))
tree.t <- as.character(sapply(names(obs),function(x) strsplit(x,split=' ')[[1]][1]))
type.t <- as.character(sapply(names(obs),function(x) strsplit(x,split=' ')[[1]][2]))
geno.t <- as.character(unlist(lapply(split(geno,paste(tree,type)),function(x) x[1])))
com.i <- do.call(rbind,lapply(obs,function(x) apply(x,2,sum)))
pbp <- tapply(pit$pb,paste(tree,type),function(x) sum(x)/length(x)) #percent pb
pbs.t <- tapply(pit$pb,tree,sum) #total across both leaf types
pbd.t <- tapply(pit$pb[type=='live'],tree[type=='live'],function(x) sum(x)/length(x)) - tapply(pit$pb[type=='sen'],tree[type=='sen'],function(x) sum(x)/length(x)) #percent difference
pbf <- tapply(pit$pb,paste(tree,type),table) #pemphigus frequencies
pbf <- do.call(rbind,lapply(pbf,function(x,n)c(x,rep(0,n-length(x))) , n=max(unlist(lapply(pbf,length)))))
pbf.liv <- do.call(rbind,lapply(lapply(split(pbf[type.t=='live',],geno.t[type.t=='live']),matrix,ncol=4),function(x) apply(x,2,sum)))
pbf.sen <- do.call(rbind,lapply(lapply(split(pbf[type.t=='sen',],geno.t[type.t=='sen']),matrix,ncol=4),function(x) apply(x,2,sum)))
rich.t <- apply(sign(com.i),1,sum)
com.rel.i <- apply(com.i,2,function(x) x/max(x))
dist.t <- as.matrix(vegdist(cbind(com.rel.i,ds=rep(min(com.rel.i[com.rel.i!=0]),nrow(com.rel.i)))))
dist.t <- dist.t[order(type.t),order(type.t)]
pdtype <- unlist(sapply(colnames(dist.t),function(x) strsplit(x,split=' ')[[1]][2]))
pd.t <- diag(dist.t[pdtype=='live',pdtype=='sen'])

## Genetic effect on P. betae
                                        #effect of leaf type
t.test(pbd.t)
cgREML(pbd.t,geno.t[type.t=='live'])
                                        #total PB across live and senescent leaves
cgREML(pbs.t,geno.t[type.t=='live'])
                                        #test of PB~genotype by type
cgREML(pbp[type.t=='live'],geno.t[type.t=='live'])
cgREML(pbp[type.t=='sen'],geno.t[type.t=='sen'])
                                        #genotype and pb frequencies on leaves
cgREML(pbf[type.t=='live',1],geno.t[type.t=='live'])
cgREML(pbf[type.t=='live',2],geno.t[type.t=='live'])
cgREML(pbf[type.t=='live',3],geno.t[type.t=='live'])
cgREML(pbf[type.t=='sen',1],geno.t[type.t=='sen'])
cgREML(pbf[type.t=='sen',2],geno.t[type.t=='sen'])
cgREML(pbf[type.t=='sen',3],geno.t[type.t=='sen'])
cgREML(pbf[type.t=='sen',4],geno.t[type.t=='sen'])
chisq.test(pbf.liv[,apply(pbf.liv,2,sum)!=0])
chisq.test(pbf.sen)

## Genotype effect on richness (i.e. individual degree)
cgREML(rich.t[type.t=='live'],geno.t[type.t=='live'])
cgREML(rich.t[type.t=='sen'],geno.t[type.t=='sen'])
cgREML((rich.t[type.t=='live']-rich.t[type.t=='sen']),geno.t[type.t=='live'])

## Genotype effect on composition
t.test(pd.t)
cgREML(pd.t,geno.t[type.t=='live'])
adonis(cbind(com.rel.i,ds=rep(min(com.rel.i[com.rel.i!=0]),nrow(com.rel.i)))~type.t)
adonis(cbind(com.rel.i,ds=rep(min(com.rel.i[com.rel.i!=0]),nrow(com.rel.i)))[type.t=='live',]~geno.t[type.t=='live'])
adonis(cbind(com.rel.i,ds=rep(min(com.rel.i[com.rel.i!=0]),nrow(com.rel.i)))[type.t=='sen',]~geno.t[type.t=='sen'])
pair.permanova(cbind(com.rel.i,ds=rep(min(com.rel.i[com.rel.i!=0]),nrow(com.rel.i)))[type.t=='live',],geno.t[type.t=='live'])$p.mat

## Genetic effect on SES
                                        #test co-occurrence across trees
                                        #values from hoth script run
                                        #SES values ~ genotype
acn.cnm <- read.csv('../data/acn_ses.csv')
###Run another ses for each tree, but put sen and live together
acn.cnm.ls <- read.csv('../data/acn_ses_ls.csv')
##
acn.ses <- acn.cnm$SES
acn.ses[is.na(acn.ses)] <- 0
cgREML(acn.ses,geno.t[order(type.t)])
cgREML(acn.ses[type.t=='live'],geno.t[type.t=='live'])
cgREML(acn.ses[type.t=='sen'],geno.t[type.t=='sen'])
cgREML(acn.cnm.ls$SES,geno.t[type.t=='live'])
#Stand level co-occurrence patterns
cnm.liv <- cnm.test(com.i[type.t=='live',])
cnm.sen <- cnm.test(com.i[type.t=='sen',])


###Building networks using tree level data
build.bpn <- function(x,alpha=0.05,p=0.05,adjust=FALSE){
  p.out <- apply(x,2,function(x) as.numeric(unlist(binom.test(sum(sign(x)),length(x),p=p))[3]))
  if (adjust){p.adjust(p.out,method='fdr')}
  x.out <- apply(sign(x),2,sum)/nrow(x)
  x.out[p.out<=alpha] <- 0
  return(x.out)
}
acn.bpn <- do.call(rbind,lapply(obs,build.bpn))
acn.type <- unlist(sapply(rownames(acn.bpn),function(x) strsplit(x,split=' ')[[1]][2]))
cgPlotweb(acn.bpn[acn.type=='live',],geno.t[type.t=='live'])
cgPlotweb(acn.bpn[acn.type=='sen',],geno.t[type.t=='live'])
                                        #tree level
                                        #run on hoth

                                        #genotype level
## gnest.liv <- list()
## gnest.liv[[1]] <- oecosimu(mean.g(acn.bpn[acn.type=='liv',],liv.geno),nestfun='nestedtemp',method='r00',alternative='greater')
## gnest.liv[[2]] <- oecosimu(mean.g(acn.bpn[acn.type=='liv',],liv.geno),nestfun='nestedtemp',method='r0',alternative='greater')
## gnest.liv[[3]] <- oecosimu(mean.g(acn.bpn[acn.type=='liv',],liv.geno),nestfun='nestedtemp',method='c0',alternative='greater')
## gnest.liv[[4]] <- oecosimu(mean.g(acn.bpn[acn.type=='liv',],liv.geno),nestfun='nestedtemp',method='r1',alternative='greater')
## gnest.sen <- list()
## gnest.sen[[1]] <- oecosimu(mean.g(acn.bpn[acn.type=='sen',],sen.geno),nestfun='nestedtemp',method='r00',alternative='greater')
## gnest.sen[[2]] <- oecosimu(mean.g(acn.bpn[acn.type=='sen',],sen.geno),nestfun='nestedtemp',method='r0',alternative='greater')
## gnest.sen[[3]] <- oecosimu(mean.g(acn.bpn[acn.type=='sen',],sen.geno),nestfun='nestedtemp',method='c0',alternative='greater')
## gnest.sen[[4]] <- oecosimu(mean.g(acn.bpn[acn.type=='sen',],sen.geno),nestfun='nestedtemp',method='r1',alternative='greater')

## Co-occurrence network structure
pit.trees <- split(pit.com,paste(leaf.type,tree))
net.trees <- lapply(pit.trees,CoNetwork)
net.d <- netDist(net.trees)
adonis(net.d~ls.geno*tree.leaf)
                                        #using the tested percent abundance values from bpn
acn.net.liv <- CoNetwork(pit.com[leaf.type=='live',])
acn.net.sen <- CoNetwork(pit.com[leaf.type=='sen',])
par(mfrow=c(1,2))
mgp(acn.net.liv,pit.com[leaf.type=='live',])
mgp(acn.net.sen,pit.com[leaf.type=='sen',])

