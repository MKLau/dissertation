###Nestedness using the bpn method to resolve bipartite network structure
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
                                        #Sourcing ComGenR on hoth
oldwd <- getwd()
setwd('~/projects/packages/ComGenR/R/')
cgn.list <- (sapply(dir(),grepl,pattern='~')|sapply(dir(),grepl,pattern='\\#'))==FALSE
sapply(dir()[cgn.list],source)
setwd(oldwd)

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
nest.liv <- list()
nest.liv[[1]] <- oecosimu(acn.bpn[acn.type=='live',],nestfun='nestedtemp',method='r00',alternative='greater',nsimul=1000)
nest.liv[[2]] <- oecosimu(acn.bpn[acn.type=='live',],nestfun='nestedtemp',method='r0',alternative='greater',nsimul=1000)
nest.liv[[3]] <- oecosimu(acn.bpn[acn.type=='live',],nestfun='nestedtemp',method='c0',alternative='greater',nsimul=1000)
nest.liv[[4]] <- oecosimu(acn.bpn[acn.type=='live',],nestfun='nestedtemp',method='r1',alternative='greater',nsimul=1000)
nest.sen <- list()
nest.sen[[1]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r00',alternative='greater',nsimul=1000)
nest.sen[[2]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r0',alternative='greater',nsimul=1000)
nest.sen[[3]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='c0',alternative='greater',nsimul=1000)
nest.sen[[4]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r1',alternative='greater',nsimul=1000)
                                        #results summary
write.csv(do.call(rbind,lapply(nest.liv,function(x)(x$'oecosimu')[c(6,2,1,3,5)])),file='../results/acn_bpnnest_live.csv')
write.csv(do.call(rbind,lapply(nest.sen,function(x)(x$'oecosimu')[c(6,2,1,3,5)])),file='../results/acn_bpnnest_sen.csv')
