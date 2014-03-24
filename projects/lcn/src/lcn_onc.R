###LCN: ONC Garden Analyses
###MKLau
###21Mar2014

library(bipartite)
library(gplots)
library(ecodist)
source('~/projects/dissertation/projects/lcn/src/lcn_load_onc.R')
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

mod.test <- oecosimu(oc,mm,method='r1',burnin=100,thin=10,nsimul=25)
mods <- computeModules(apply(oc[,-1],2,function(x) x/max(x)))
plotModuleWeb(mods)
                                        #bipartite network
plotweb(oc[order(apply(oc,1,sum),decreasing=TRUE),order(apply(oc,2,sum),decreasing=TRUE)],text.rot=90,method='normal')
plotweb(omu[order(apply(omu,1,sum),decreasing=TRUE),order(apply(omu,2,sum),decreasing=TRUE)],
        text.rot=90,method='normal',labsize=2)

## SES ~ Roughness
summary(lm(sqrt(abs(os))~prb))

## SES ~ genotype
plot(os~prb,xlab='Percent Rough Bark',ylab='SES',pch=19)
abline((lm(os~prb)))
barplot2(osgmu[order(osgmu,decreasing=TRUE)],plot.ci=TRUE,
         ci.u=osgmu[order(osgmu,decreasing=TRUE)]+osgse[order(osgmu,decreasing=TRUE)],
         ci.l=osgmu[order(osgmu,decreasing=TRUE)]-osgse[order(osgmu,decreasing=TRUE)],
         ylab='SES')
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

## co ~ Roughness
## Specific co-occurrence pairs are not changing with bark roughness
adonis(cbind(oco,ds=rep(1,nrow(oco)))~prb,permutations=5000)
adonis(cbind(apply(oco,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(oco)))~prb,permutations=5000)

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

## Araujo network
oan <- coNet(do.call(rbind,oq))
oan <- log(oan+1)
oan.col <- coNet(do.call(rbind,oq),showsign=TRUE)
oan.vs <- (apply(oc,2,sum))
oan.vs <- oan.vs/max(oan.vs)+1.5
gplot(oan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='grey',
      vertex.cex=oan.vs,edge.lwd=oan)
text(coord,rownames(oan))

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
