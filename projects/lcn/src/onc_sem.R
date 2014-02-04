###LCO - SEM for lichen at ONC
###Taken out of the notebook.Rnw file chunk
###4 Feb 2014

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
                                        #SEM
library(sem)
sem.geno <- array(rep(onc.geno,length(unique(onc.geno))),dim=c(length(onc.geno),length(unique(onc.geno))))
for (i in 1:length(unique(onc.geno))){
  sem.geno[sem.geno[,i]!=unique(onc.geno)[i],i] <- 0
  sem.geno[sem.geno[,i]==unique(onc.geno)[i],i] <- 1
}
sem.geno <- apply(sem.geno,2,as.numeric)
colnames(sem.geno) <- paste('Geno',1:length(unique(onc.geno)),sep='')
sem.geno <- sem.geno[,-ncol(sem.geno)]
sem.data <- cbind(sem.geno,roughness=onc.rough,Axis1=onc.rot)
model.sem <- specifyModel(file='../../lcn/data/sem/onc_sem.txt')
                                        #transforms???
                                        #sem fitting and analysis
                                        #colnames(sem.data)
Sigma <- cov(sem.data)
sem.fit <- sem(model.sem,S=Sigma,N=nrow(sem.data))
summary(sem.fit)
modIndices(sem.fit)
                                        #effects(sem.fit)
                                        #hist(residuals(sem.fit))
                                        #stdCoef(sem.fit)
pathDiagram(sem.fit,file='../../lcn/results/onc_sem.png',edge.labels='values',standardize=TRUE
            ,ignore.double=TRUE,size=c(12,12),edge.font=c("Arial", 10),
            graphics.fmt='png') #export to graphviz
                                        #
