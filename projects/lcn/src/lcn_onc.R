###LCN: ONC Garden Analyses
###MKLau
###21Mar2014

rm(list=ls())
library(bipartite)
library(gplots)
library(ecodist)
cbp <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
         "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
source('~/projects/dissertation/projects/lcn/src/lcn_load_onc.R')
source('/Users/Aeolus/projects/dissertation/projects/lcn/src/getMods.R')
source('~/projects/packages/cooc/src/cooc.R')
source('~/projects/packages/ComGenR/R/coNet.R')

###Data objects:
###oc = co-occurrences summed across all cells for each tree
###oq = co-occurrence matrices separated out for each tree
###og = genotypes
###osgmu = ses genotype means
###osgse = ses genotype SE
###os = onc ses values
###oco = co-occurrence counts
###och = checker counts
###prb = percent rough bark (averaged between the upper and lower)
###Data notes:
###No physciods
###Lecanoras merged

## Height composition
### Previous tests show that this is not significant

## Overall Co-occurrence Patterns
#oecosimu(oc,cs,method='r1',burnin=100,thin=10,nsimul=5000)
## Call: oecosimu(comm = oc, nestfun = cs, method = "r1", nsimul = 5000,burnin = 100, thin = 10)
## simulation method r1 with 5000 simulations
##           statistic       z    mean    2.5%     50%  97.5% Pr(sim.)    
## statistic    26.333 -6.3127  88.948  69.972  88.944 108.81    2e-04

## Nestedness, Modularity, Centrality
## oecosimu(oc,nested,method='r1',burnin=100,thin=10,nsimul=5000)
##             statistic       z   mean   2.5%    50%  97.5% Pr(sim.)    
## binmatnest2    8.4664 -4.8652 23.999 17.912 23.902 30.427    2e-04
## mod.test <- oecosimu(apply(omu,2,function(x) x/max(x)),mm,method='r1',burnin=100,thin=10,nsimul=25)
## mods <- computeModules(apply(omu,2,function(x) x/max(x)))
## oom <- slot(mods,'likelihood')
## omp <- length(mod.test$oecosimu$simulated[mod.test$oecosimu$simulated>=oom])/length(mod.test$oecosimu$simulated)
## omz <- (oom - mean(mod.test$oecosimu$simulated)) / sd(mod.test$oecosimu$simulated)
### obs = 0.20,z=9.637778,p<0.001
                                        #bipartite network
## mod.col <- getMods(mods)
## col.trees <- mod.col[[1]]
## col.spp <- mod.col[[2]]
## col.trees[mod.col[[1]]==3] <- 2
## col.spp[mod.col[[2]]==3] <- 2
## col.trees[mod.col[[1]]==2] <- 3
## col.spp[mod.col[[2]]==2] <- 3
## col.trees <- cbp[col.trees]
## col.spp <- cbp[col.spp]
## plotweb(omu[order(apply(omu,1,sum),decreasing=TRUE),order(apply(omu,2,sum),decreasing=TRUE)],
##         text.rot=90,method='normal',labsize=2,
##         col.low=col.trees[order(apply(omu,1,sum),decreasing=TRUE)],
##         col.high=col.spp[order(apply(omu,2,sum),decreasing=TRUE)])

## SES ~ Roughness
summary(rlm(sqrt(abs(os))~prb))

###transform percent rough bark too
cor.test(sqrt(abs(os)),sqrt(prb))
##t = 1.2843, df = 55, p-value = 0.2044

###Mantel
os.d <- dist(as.matrix(os))
prb.d <- dist(as.matrix(prb))
mantel(os.d~prb.d)

##    mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
##-0.08334190  0.97800000  0.02300000  0.06900000 -0.10906468 -0.06058215 
os.d <- dist(sqrt(abs(as.matrix(tapply(os,og,mean)))))
prb.d <- dist(as.matrix(tapply(prb,og,mean)))
mantel(os.d~prb.d)


mantel os.d prb.d
library(sem)
###OS <- PRB <- OG
##make genotype a factor
og.mat <- array(0,dim=c(length(og),length(unique(og))))
colnames(og.mat) <- paste('Geno',1:ncol(og.mat),sep='')
for (i in 1:length(unique(og))){
  og. <- og
  og.[og==unique(og)[i]] <- 1
  og.[og!=unique(og)[i]] <- 0
  og.mat[,i] <- as.numeric(og.)
}

library(sem)
   model.sem <- specifyModel(file='../data/sem_os.txt')
   sem.data <- data.frame(og.mat[,unique(og)!='1005'],prb,os)
                                        #sem.data[,1:12] <- apply(sem.data[,1:12],2,factor)
                                        #transforms
   os <- sqrt(abs(os))
                                        #sem fitting and analysis
                                        #colnames(sem.data)
   Sigma <- var(sem.data)
   sem.fit <- sem(model.sem,S=Sigma,N=nrow(sem.data))
   summary(sem.fit)
   modIndices(sem.fit)
                                        #effects(sem.fit)
                                        #hist(residuals(sem.fit))
                                        #stdCoef(sem.fit)
   pathDiagram(sem.fit,file='semPathA',edge.labels='values',standardize=TRUE
               ,ignore.double=FALSE,size=c(12,12),edge.font=c("Arial", 10),
               graphics.fmt='png') #export to graphviz



## SES ~ genotype
oneway.test(os~og)
plot(os~prb,xlab='Percent Rough Bark',ylab='SES',pch=19)
abline((lm(os~prb)))
barplot2(osgmu[order(osgmu,decreasing=TRUE)],plot.ci=TRUE,
         ci.u=osgmu[order(osgmu,decreasing=TRUE)]+osgse[order(osgmu,decreasing=TRUE)],
         ci.l=osgmu[order(osgmu,decreasing=TRUE)]-osgse[order(osgmu,decreasing=TRUE)],
         ylab='SES')

### PRB~genotype
oneway.test(prb~og)
prb.mu <- tapply(prb,og,mean)
prb.se <- tapply(prb,og,function(x) sd(x)/length(x))
barplot2(prb.mu[order(prb.mu,decreasing=TRUE)],plot.ci=TRUE,
         ci.u=prb.mu[order(prb.mu,decreasing=TRUE)]+prb.se[order(prb.mu,decreasing=TRUE)],
         ci.l=prb.mu[order(prb.mu,decreasing=TRUE)]-prb.se[order(prb.mu,decreasing=TRUE)],
         ylab='Percent Rough Bark')
                                        #beside plot
mu <- rbind('Percent Rough Bark'=prb.mu,'SES'=osgmu[match(names(prb.mu),names(osgmu))])
se <- rbind('prb.se'=prb.se,'ses.se'=osgse[match(names(prb.se),names(osgse))])
barplot2(mu[,order(mu[2,],decreasing=TRUE)],plot.ci=TRUE,
         ci.u=mu[,order(mu[2,],decreasing=TRUE)]+se[,order(mu[2,],decreasing=TRUE)],
         ci.l=mu[,order(mu[2,],decreasing=TRUE)]-se[,order(mu[2,],decreasing=TRUE)],
         ylab='',beside=TRUE)
                                        #cross-hair plot
ch.x <- cbind('Percent Rough Bark'=prb,'SES'=os)
ch.plot(ch.x,og)
abline(lm(osgmu~prb.mu),lty=2)
                                        #ses ~ rflp distance
mantel(oms.d~rflp.d)
mantel(oms.d~dist(oprbmu))
mantel(dist(oprbmu)~rflp.d)
plot(oms.d~rflp.d,pch=19,xlab='RFLP Genetic Relatedness',ylab='SES Dissimilarity')
abline(lm(oms.d~rflp.d))

plot(oms.d~dist(oprbmu),pch=19,xlab='Bark Roughness Dissimilarity',ylab='SES Dissimilarity')
abline(lm(oms.d~dist(oprbmu)))
mantel(oms.d~dist(oprbmu))
                                        #
y <- (oms[match(rownames(as.matrix(rflp.d)),names(oms))])
x <- oprbmu[match(rownames(as.matrix(rflp.d)),names(oprbmu))]
summary(lm(y~x))
plot(y~x,xlab='Mean Bark Roughness',ylab='Mean SES',pch=19)
abline(lm(y~x))


## Araujo network
oan <- coNet(do.call(rbind,oq))
oan <- log(oan+1)
oan.col <- coNet(do.call(rbind,oq),showsign=TRUE)
oan.vs <- (apply(oc,2,sum))
oan.vs <- oan.vs/max(oan.vs)+1.5
gplot(oan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='lightgrey',
      vertex.cex=oan.vs,edge.lwd=oan)
text(coord,rownames(oan))

## co ~ Roughness
## Specific co-occurrence pairs are not changing with bark roughness
sigCo <- function(x){
  y <- character()
  k <- 0
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (i!=j&(x[i,j]!=0)){
        k <- k + 1
        y[k] <- paste(colnames(x)[i],colnames(x)[j],sep='_')
      }
    }
  }
  return(y)
}
                                        #
adonis(cbind(oco,ds=rep(1,nrow(oco)))~prb,permutations=5000)
adonis(cbind(apply(oco,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(oco)))~prb,permutations=5000)
oco.sig <- cbind(oco[,colnames(oco)%in%sigCo(oan)],ds=rep(1,nrow(oco)))
och.sig <- cbind(och[,colnames(och)%in%sigCo(oan)],ds=rep(1,nrow(och)))
oco.sig.mu <- apply(oco.sig,2,function(x,g) tapply(x,g,mean),g=og)
och.sig.mu <- apply(och.sig,2,function(x,g) tapply(x,g,mean),g=og)
oco.sig.mu.d <- as.matrix(vegdist(oco.sig.mu))
oco.sig.mu.d <- oco.sig.mu.d[match(rownames(as.matrix(rflp.d)),rownames(oco.sig.mu.d)),
                             match(rownames(as.matrix(rflp.d)),rownames(oco.sig.mu.d))]
oco.sig.mu.d <- as.dist(oco.sig.mu.d)
och.sig.mu.d <- as.matrix(vegdist(och.sig.mu))
och.sig.mu.d <- och.sig.mu.d[match(rownames(as.matrix(rflp.d)),rownames(och.sig.mu.d)),
                             match(rownames(as.matrix(rflp.d)),rownames(och.sig.mu.d))]
och.sig.mu.d <- as.dist(och.sig.mu.d)
                                        #
adonis(oco.sig~prb)
adonis(oco.sig~og)
adonis(och.sig~prb)
adonis(och.sig~og)
mantel(oco.sig.mu.d~rflp.d)
mantel(och.sig.mu.d~rflp.d)

##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## prb        1    0.5294 0.52941  1.4714 0.02606  0.156
## Residuals 55   19.7894 0.35981         0.97394       
## Total     56   20.3188                 1.00000       
adonis(cbind(oco,ds=rep(1,nrow(oco)))~og,permutations=5000)
adonis(cbind(apply(oco,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(oco)))~og,permutations=5000)

## ch ~ Roughness
## Pairs of checkers are changing with roughness
adonis(cbind(och,ds=rep(1,nrow(och)))~prb,permutations=5000)
adonis(cbind(apply(och,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(och)))~prb,permutations=5000)
adonis(cbind(och,ds=rep(1,nrow(och)))~og,permutations=5000)
adonis(cbind(apply(och,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(och)))~og,permutations=5000)

## ch ~ co ~ SES
nco <- apply(oco,1,sum)
nch <- apply(och,1,sum)
adonis(cbind(oco,ds=rep(1,nrow(oco)))~os,permutations=5000)
adonis(cbind(och,ds=rep(1,nrow(och)))~os,permutations=5000)
mantel(vegdist(cbind(oco,ds=rep(0.0001,nrow(oco)))),vegdist(och))
summary(lm(apply(oco,2,sum)~apply(och,2,sum)))
summary(lm(nco~nch))
                                        #ch~ses
plot(dist(as.matrix(os))~vegdist(och),xlab='Checker Dissimilarity',
     ylab='SES Dissimilarity',pch=19)
abline(lm(dist(as.matrix(os))~vegdist(och)))
                                        #
plot(vegdist(cbind(oco,ds=rep(1,nrow(oco))))~vegdist(cbind(och,ds=rep(1,nrow(och)))),pch=19,
     xlab='Checker Dissimilarity',ylab='Co-Occurrence Dissimilarity')
abline(lm(vegdist(cbind(oco,ds=rep(1,nrow(oco))))~vegdist(cbind(och,ds=rep(1,nrow(och))))))


### RFLP distance ordination with overlay
rflp.ord <- princomp(rflp.d)
rflp.ord$scores[,1:2] <- rflp.ord$scores[,2:1]
vector.in <- cbind(oms[match(rownames(as.matrix(rflp.d)),names(oms))],oprbmu)
colnames(vector.in) <- c('SES','PRB')
rflp.vec <- envfit(rflp.ord$scores[,1:2],vector.in)
                                        #jitter 1017
rflp.ord$scores[rownames(rflp.ord$scores)=='1017',1:2] <- rflp.ord$scores[rownames(rflp.ord$scores)=='1017',1:2] + 0.05

plot(rflp.ord$scores[,1:2],pch='',xlab='PCA 1',ylab='PCA 2')
text(rflp.ord$scores[,1:2],labels=rownames(as.matrix(rflp.d)))
plot(rflp.vec,col=2)


###average network
oaan <- lapply(oq,coNet)
oaan <- Reduce('+',oaan)
oaan <- log(oaan+1)
sum(abs(sign(oaan)-sign(oan)))

gplot(oan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='lightgrey',
      vertex.cex=oan.vs,edge.lwd=oan)
text(coord,rownames(oan))

gplot(oaan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='lightgrey',
      vertex.cex=oan.vs,edge.lwd=oaan)
text(coord,rownames(oaan))

###Eigen Centrality
library(igraph)
onec <- centralization.evcent(graph.adjacency(oan,weighted=TRUE))$vector
plot(graph.adjacency(oan,weighted=TRUE),vertex.size=onec*10)
names(onec) <- rownames(oan)
