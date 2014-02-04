###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

rm(list=ls())
library(ComGenR)
library(lme4)
cgREML <- function(x,g,cv='covariate'){
  if (any(cv=='covariate')){
    x.lmer <- lmer(x~(1|g))
    x.lm <- lm(x~1)
    chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
    p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  }else{
    x.lmer <- lmer(x~(1|g)+cv)
    x.lm <- lm(x~1)
    chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
    p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  }
  return(c(chi2=chi2,P.value=p.chi2))
}
build.bpn <- function(x,alpha=0.05,p=0.001,adjust=FALSE){
  p.out <- apply(x,2,function(x) as.numeric(unlist(binom.test(sum(sign(x)),length(x),p=p))[3]))
  if (adjust){p.adjust(p.out,method='fdr')}
  x.out <- apply(sign(x),2,sum)/nrow(x)
  x.out[p.out>alpha] <- 0
  return(x.out)
}

###Garden Analysis
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
obs.cs <- unlist(lapply(onc.q,cscore))
                                        #load ses values
onc.ses <- dget('../../lcn/data/onc_tree_ses.Rdata')
onc.p <- dget('../../lcn/data/onc_tree_pval.Rdata')
onc.tn <- lapply(onc.q,CoNetwork) #tree level networks
if (any(names(onc.q)=='')){names(onc.ses) <- names(onc.q)[names(onc.q)!='']}
onc.ses[is.na(onc.ses)] <- 0
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
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
onc.rough <- avg.rough[match(ses.tree,r.tree)]
all(ses.tree==names(onc.rough))
                                        #Microsat data from Nash
gen.d <- read.csv(file='../../lcn/data/ONC_MSAT_datafromnash.csv')[,-1]
gen.d[is.na(gen.d)] <- 0
gen.d <- as.dist((gen.d))
                                        #community data
onc.com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
onc.R <- apply(sign(onc.com),1,sum)
onc.H <- diversity(onc.com)
onc.com.rel <- apply(onc.com,2,function(x) x/max(x))
onc.com.rel <- cbind(onc.com.rel,ds=rep(min(onc.com.rel[onc.com.rel!=0]),nrow(onc.com.rel)))
onc.com <- cbind(onc.com,ds=rep(min(onc.com[onc.com!=0]),nrow(onc.com)))
                                        #nmds and procrustes rotation
onc.nms <- nmds.min(nmds(vegdist(onc.com.rel),2,2))
onc.rot <- procrustes(dist(onc.rough),onc.nms,scale=TRUE)
onc.rot <- t(onc.rot$rotation)
onc.rot.sem <- onc.rot
onc.rot <- onc.rot[,(1:2)[abs(cor(cbind(onc.rough,onc.rot))[1,2:3])==max(abs(cor(cbind(onc.rough,onc.rot))[1,2:3]))]]

###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env);adonis(com~height);detach(env)

##Species accumulation curve
sac.com <- list()
for (i in 1:length(unique(env$tree))){
  sac.com[[i]] <- apply(com[env$tree==unique(env$tree)[i],],2,sum)
}
sac.com <- do.call(rbind,sac.com)
plot(specaccum(sac.com),ylim=c(0,12),ylab='Species Richness',xlab='Trees Sampled',font.lab=2)
legend('topright',legend=c('Garden','Wild'),lty=c(1,1),col=c(1,2))

###Roughness effects on composition and richness
summary(lm(onc.R~onc.rough))
summary(lm(onc.H~onc.rough))
summary(lm(onc.rot~onc.rough))
adonis(onc.com~onc.rough,permutations=5000)
adonis(onc.com.rel~onc.rough,permutations=5000)

###Genotype effects on roughness, composition, richness
cgREML(onc.rough,onc.geno)
getH2C(onc.rough,onc.geno,method='nmds')
cgREML(onc.R,onc.geno)
cgREML(onc.H,onc.geno)
cgREML(onc.rot,onc.geno)

###modeling
###Co-occurrences
                                        #co-occurrence
stand.null <- unlist(dget(file='../../lcn/data/onc_stand_null.Rdata'))
stand.ses <- (cscore(onc[,7:ncol(onc)]) - mean(stand.null)) / sd(stand.null)
stand.ses.p <- length(stand.null[stand.null<=stand.ses])/length(stand.null)
c(stand.ses,stand.ses.p)
                                        #
onc.treemean.ses <- dget(file='../../lcn/results/onc_ses_treemean.Rdata')
onc.treemean.ses
                                        #
onc.ses <- dget(file='../../lcn/data/onc_tree_ses.Rdata')
onc.ses[is.na(onc.ses)] <- 0
cgREML(onc.ses,onc.geno)
plot(onc.ses~as.numeric(factor(onc.geno)))
summary(lm(onc.ses~onc.rough))
onc.mg.ses <- tapply(onc.ses,onc.geno,mean)
onc.mg.rough <- tapply(onc.rough,onc.geno,mean)
plot(onc.mg.ses^2~onc.mg.rough)
abline(lm(onc.mg.ses^2~onc.mg.rough))
summary(lm((onc.mg.ses^2)~onc.mg.rough))
                                        #unipartite networks
lcn <- CoNetwork((onc[,7:ncol(onc)]))
mgp(lcn,(onc[,7:ncol(onc)]),displaylabels=TRUE,loc=FALSE)
                                        #bipartite network
onc.bpn <- do.call(rbind,lapply(onc.q,build.bpn))
cgPlotweb(onc.bpn,onc.geno)
## onc.nest <- list()
## onc.nest[[1]] <- oecosimu(onc.bpn,nestfun='nestedtemp',method='r00',nsimul=1000)
## onc.nest[[2]] <- oecosimu(onc.bpn,nestfun='nestedtemp',method='r0',nsimul=1000)
## onc.nest[[3]] <- oecosimu(onc.bpn,nestfun='nestedtemp',method='c0',nsimul=1000)
## onc.nest[[4]] <- oecosimu(onc.bpn,nestfun='nestedtemp',method='r1',nsimul=1000)
## onc.geno.nest <- list()
## onc.geno.nest[[1]] <- oecosimu(mean.g(onc.bpn,onc.geno),nestfun='nestedtemp',method='r00',nsimul=1000)
## onc.geno.nest[[2]] <- oecosimu(mean.g(onc.bpn,onc.geno),nestfun='nestedtemp',method='r0',nsimul=1000)
## onc.geno.nest[[3]] <- oecosimu(mean.g(onc.bpn,onc.geno),nestfun='nestedtemp',method='c0',nsimul=1000)
## onc.geno.nest[[4]] <- oecosimu(mean.g(onc.bpn,onc.geno),nestfun='nestedtemp',method='r1',nsimul=1000)
onc.nest <- read.csv(file='../../lcn/results/onc_nest_tree.csv')
onc.geno.nest <- read.csv(file='../../lcn/results/onc_nest_geno.csv')
                                        #correlate genotype mean roughness and mean degree
onc.deg <- apply(sign((onc.bpn)),1,sum)
summary(lm(onc.deg~onc.rough))
