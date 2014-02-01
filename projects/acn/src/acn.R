###Final set of analyses for the ACN paper
###29 Jan 2014
###To be run up on Hoth

###Re-do labeling so that any thing at the:
##cell scale is labeled c
##quadrat scale is labeled q
##tree scale is labeled t

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
tree <- pit$tree
liv.geno <- as.character(unlist(lapply(split(geno$liv,split(tree,pit[,1])$liv),function(x) x[1])))
sen.geno <- as.character(unlist(lapply(split(geno$sen,split(tree,pit[,1])$sen),function(x) x[1])))
liv.com <- do.call(rbind,lapply(split(liv,split(tree,pit[,1])$liv),function(x) apply(x,2,sum)))
sen.com <- do.call(rbind,lapply(split(sen,split(tree,pit[,1])$sen),function(x) apply(x,2,sum)))
tree.geno <- as.character(unlist(lapply(split(unlist(geno),tree),function(x) x[1])))
acn.geno <- c(liv.geno,sen.geno)
leaf.type <- as.character(pit[,1])
tree.leaf <- unlist(lapply(split(leaf.type,paste(leaf.type,tree)),function(x) x[1]))
n.leaves <- unlist(lapply(split(pit.com,paste(leaf.type,tree)),function(x) nrow(x)))
liv.pc <- liv.com
for (i in 1:nrow(liv.com)){liv.pc[i,] <- liv.com[i,]/split(n.leaves,tree.leaf)[[1]][i]}
sen.pc <- sen.com
for (i in 1:nrow(sen.com)){sen.pc[i,] <- sen.com[i,]/split(n.leaves,tree.leaf)[[2]][i]}
pit.trees <- split(pit.com,paste(tree,leaf.type))

## Genetic effect on P. betae
                                        #total PB across live and senescent leaves
pb.total <- tapply(pit$pb,tree,sum)
summary(aov(pb.total~tree.geno))
                                        #reml
                                        #using notes from the following link
                                        #http://www.stat.wisc.edu/~ane/st572/notes/lec21.pdf
library(lme4)
cgREML <- function(x,g){
  x.lmer <- lmer(x~(1|g))
  x.lm <- lm(x~1)
  chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
  p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  return(c(chi2=chi2,P.value=p.chi2))
}
ls.pbd <- liv.pc[,colnames(liv.pc)=='pb']-sen.pc[,colnames(sen.pc)=='pb']
                                        #test of PB as percent leaves serveyed
cgREML(pb.total,tree.geno)
cgREML(liv.com[,colnames(liv.com)=='pb'],tree.geno)
cgREML(sen.com[,colnames(sen.com)=='pb'],tree.geno)
cgREML(liv.pc[,colnames(liv.pc)=='pb'],tree.geno)
cgREML(sen.pc[,colnames(sen.pc)=='pb'],tree.geno)
t.test(ls.pbd)
cgREML(ls.pbd,tree.geno)


## Genotype effect on composition
                                        #paired test
ls.com <- rbind(liv.com,sen.com)
ls.com.rel <- apply(ls.com,2,function(x) if(all(x==0)){x}else{x/max(x)})
ls.com <- cbind(ls.com,ds=rep(min(ls.com[ls.com!=0],nrow(ls.com))))
ls.com.rel <- cbind(ls.com.rel,ds=rep(min(ls.com.rel[ls.com.rel!=0],nrow(ls.com.rel))))
ls.d <- vegdist(ls.com)
ls.d.rel <- vegdist(ls.com.rel)
ls.pair <- diag(as.matrix(ls.d)[tree.leaf=='live',tree.leaf=='sen'])
ls.pair.rel <- diag(as.matrix(ls.d.rel)[tree.leaf=='live',tree.leaf=='sen'])
liv.com.rel <- apply(liv.com,2,function(x) x/max(x))
liv.com.rel[is.na(liv.com.rel)] <- 0
sen.com. <- cbind(sen.com,ds=rep(1,nrow(sen.com)))
sen.com.rel <- apply(sen.com.,2,function(x) x/max(x))
sen.com.rel[,ncol(sen.com.rel)] <- min(sen.com.rel[sen.com.rel!=0])
liv.pc.rel <- apply(liv.pc,2,function(x) if(all(x==0)){x}else{x/max(x)})
sen.pc.rel <- apply(sen.pc,2,function(x) if(all(x==0)){x}else{x/max(x)})
liv.pc.rel <- cbind(liv.pc.rel,ds=rep(min(liv.pc.rel[liv.pc.rel!=0]),nrow(liv.pc.rel)))
sen.pc.rel <- cbind(sen.pc.rel,ds=rep(min(sen.pc.rel[sen.pc.rel!=0]),nrow(sen.pc.rel)))
ls.pc <- rbind(liv.pc,sen.pc)
ls.pc <- cbind(ls.pc,ds=rep(min(ls.pc[ls.pc!=0]),nrow(ls.pc)))
ls.geno <- c(liv.geno,sen.geno)
                                        #paired tests
t.test(ls.pair)
t.test(ls.pair.rel)
cgREML(ls.pair,tree.geno)
cgREML(ls.pair.rel,tree.geno)
                                        #non-paired test
adonis(liv.com~liv.geno)
adonis(sen.com.~sen.geno)
adonis(liv.com.rel~liv.geno)
adonis(sen.com.rel~sen.geno)
adonis(ls.pc~ls.geno*tree.leaf)
adonis(ls.pc.rel~ls.geno*tree.leaf)
adonis(liv.pc.rel~liv.geno)
adonis(sen.pc.rel~sen.geno)

## Genetic effect on SES
all.ses <- cnm.test(rbind(liv.com,sen.com),nits=1000)
test <- list()
for (i in 1:10){test[[i]] <- cnm.test(rbind(liv.com[sample(1:nrow(liv.com),17),],sen.com[sample(1:nrow(liv.com),17),]),nits=1000)}
liv.ses <- cnm.test(liv.com,nits=1000)
sen.ses <- cnm.test(sen.com,nits=1000)
rbind(all.ses,test.ses=apply(do.call(rbind,test),2,mean),liv.ses,sen.ses)
                                        #test co-occurrence across trees
                                        #values from hoth script run
                                        #SES values ~ genotype
acn.cnm <- read.csv('../data/acn_ses.csv')
###Run another ses for each tree, but put sen and live together
###call it acn_ses_ls.csv
acn.ses <- acn.cnm$SES
acn.ses[is.na(acn.ses)] <- 0
cgREML(acn.ses,acn.geno)
cgREML(acn.ses[tree.leaf=='live'],acn.geno[tree.leaf=='live'])
cgREML(acn.ses[tree.leaf=='sen'],acn.geno[tree.leaf=='sen'])
summary(aov(acn.ses~acn.geno*tree.leaf))

###Building networks using tree level data
build.bpn <- function(x,alpha=0.05,p=0.05,adjust=FALSE){
  p.out <- apply(x,2,function(x) as.numeric(unlist(binom.test(sum(sign(x)),length(x),p=p))[3]))
  if (adjust){p.adjust(p.out,method='fdr')}
  x.out <- apply(sign(x),2,sum)/nrow(x)
  x.out[p.out<=alpha] <- 0
  return(x.out)
}
acn.bpn <- do.call(rbind,lapply(pit.trees,build.bpn))
acn.type <- substr(names(pit.trees),1,3)
cgPlotweb(acn.bpn[acn.type=='liv',],liv.geno)
cgPlotweb(acn.bpn[acn.type=='sen',],sen.geno)
                                        #tree level
nest.liv <- list()
nest.liv[[1]] <- oecosimu(acn.bpn[acn.type=='liv',],nestfun='nestedtemp',method='r00',alternative='greater',nsimul=1000)
nest.liv[[2]] <- oecosimu(acn.bpn[acn.type=='liv',],nestfun='nestedtemp',method='r0',alternative='greater',nsimul=1000)
nest.liv[[3]] <- oecosimu(acn.bpn[acn.type=='liv',],nestfun='nestedtemp',method='c0',alternative='greater',nsimul=1000)
nest.liv[[4]] <- oecosimu(acn.pbn[acn.type=='liv',],nestfun='nestedtemp',method='r1',alternative='greater',nsimul=1000)
nest.sen <- list()
nest.sen[[1]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r00',alternative='greater',nsimul=1000)
nest.sen[[2]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r0',alternative='greater',nsimul=1000)
nest.sen[[3]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='c0',alternative='greater',nsimul=1000)
nest.sen[[4]] <- oecosimu(acn.bpn[acn.type=='sen',],nestfun='nestedtemp',method='r1',alternative='greater',nsimul=1000)
                                        #results summary
do.call(rbind,lapply(nest.liv,function(x)(x$'oecosimu')[c(6,2,1,3,5)]))
do.call(rbind,lapply(nest.sen,function(x)(x$'oecosimu')[c(6,2,1,3,5)]))
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

