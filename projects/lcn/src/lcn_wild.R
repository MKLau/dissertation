###LCN: Wild Stand Analyses
###MKLau
###21Mar2014

###Meta
###Site = Uintah, UT
###Study area = 225 * 463 = 104,175 m2 = 0.104175 km2
rm(list=ls())
cbp <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
         "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
source('/Users/Aeolus/projects/dissertation/projects/lcn/docs/LCO_analyses/source/pairs.R')
source('~/projects/dissertation/projects/lcn/src/lcn_load_wild.R')
source('~/projects/packages/cooc/src/cooc.R')

###Data objects:
###wc = co-occurrences summed across all cells for each tree
###wq = co-occurrence matrices separated out for each tree
###ws = matrix of ses and related stats for trees
###wco = co-occurrence counts
###wch = checker counts
###prb = percent rough bark (averaged between the upper and lower)
###Data notes:
###No physciods
###Lecanoras merged
###Folios grey black was deleted as it occurred in only three cells on one tree

## Height composition
### Previous tests show that this is not significant

## Overall Co-occurrence Patterns
#wco.all <- oecosimu(do.call(rbind,wq),cs,method='r1',burnin=100,thin=10,nsimul=5000)
## oecosimu(wc,cs,method='r1',burnin=100,thin=10,nsimul=5000)
## statistic    1.0833 -3.1515 4.7057 2.4722 4.7222 6.9444   0.0018 **

## Nestedness, Modularity, Centrality
## oecosimu(wc,nested,method='r1',burnin=100,thin=10,nsimul=5000)
##             statistic       z   mean   2.5%    50%  97.5%  Pr(sim.)    
## binmatnest2    6.7749 -2.7125 18.127 10.389 18.018 26.648 0.0005999 ***

##wmods <- oecosimu(apply(wc[,colnames(wc)!='Pu'],2,function(x) x/max(x)),mm,method='r1',burnin=100,thin=10,nsimul=20)
## wmp <- length(wmods$oecosimu$simulated[wmods$oecosimu$simulated>=slot(wom,'likelihood')])/length(wmods$oecosimu$simulated)
## wmz <- (slot(wom,'likelihood') - mean(wmods$oecosimu$simulated))/sd(wmods$oecosimu$simulated)
##       statistic  z        mean         Pr(sim)
## mm    0.3243    7.422347   0.2099011  0.000000000
## plotModuleWeb(wom)
#wom <- computeModules(apply(wc[,colnames(wc)!='Pu'],2,function(x) x/max(x)))
wom. <- getMods(wom)
col.trees <- cbp[wom.[[1]]]
col.spp <- cbp[c(wom.[[2]][1:5],Pu=5,wom.[[2]][6:8])]
col.spp[6] <- 'grey'
rownames(wc) <- paste('Unitah',1:nrow(wc),sep='')
plotweb(wc[order(apply(wc,1,sum),decreasing=TRUE),order(apply(wc,2,sum),decreasing=TRUE)],
        text.rot=90,method='normal',labsize=1.65,
        col.low=col.trees[order(apply(wc,1,sum),decreasing=TRUE)],
        col.high=col.spp[order(apply(wc,2,sum),decreasing=TRUE)])
## SES ~ Roughness
library(MASS)
summary(rlm(ws$wses~prb))
## rlm = t=-2.1766, p=0.025
##Mantel for ses and roughness
ws.d <- dist(as.matrix(ws$wses))
prb.d <- dist(as.matrix(prb))
mantel(ws.d~prb.d)

## PRB ~ age
## summary(rlm(prb~age))
## rlm = t=2.1461, p=0.02650551
## SES ~ age
summary(rlm(ws$wses~age))
## rlm = t=-2.1766,p=

##SES ~ prb * age
summary(lm(ws$wses~prb))

plot(ws$wses~prb,xlab='Percent Rough Bark',ylab='SES',pch=19)
abline((lm(ws$wses~prb)),lty=2)

## co ~ Roughness
## Specific co-occurrence pairs are not changing with bark roughness
## adonis(cbind(apply(wco,2,function(x) if (all(x==0)){x}else{x/max(x)}),ds=rep(0.001,nrow(wco)))~prb,permutations=5000)
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## prb        1    0.2249 0.22493 0.75984 0.05955 0.5975
## Residuals 12    3.5522 0.29602         0.94045       
## Total     13    3.7771                 1.00000       

## ch ~ Roughness
## Pairs of checkers are changing with roughness
## adonis(apply(wch,2,function(x) if(all(x==0)){x}else{x/max(x)})~prb,permutations=5000)
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## prb        1   0.35875 0.35875  1.9225 0.13809 0.1058
## Residuals 12   2.23925 0.18660         0.86191
## Total     13   2.59801                 1.00000

## ch ~ co ~ SES
nco <- apply(wco,1,sum)
nch <- apply(wch,1,sum)
adonis(cbind(wco,ds=rep(1,nrow(wco)))~ws$wses,permutations=5000)
adonis(wch~ws$wses,permutations=5000)
mantel(vegdist(cbind(wco,ds=rep(0.0001,nrow(wco))))~vegdist(wch))
summary(lm(apply(wco,2,sum)~apply(wch,2,sum)))
                                        #co~ch
plot(vegdist(cbind(wco,ds=rep(0.0001,nrow(wco))))~vegdist(wch),ylab='Co-Occurrence Dissimilarity',
     xlab='Checker Dissimilarity',pch=19)
abline(lm((vegdist(cbind(wco,ds=rep(0.0001,nrow(wco))))~vegdist(wch))))
plot(apply(wco,2,sum)~apply(wch,2,sum))
abline(lm(apply(wco,2,sum)~apply(wch,2,sum)))
                                        #co~ses
plot(dist(as.matrix(ws$wses))~vegdist(cbind(wco,ds=rep(0.0001,nrow(wco)))),xlab='Co-Occurrence Dissimilarity',
     ylab='SES Dissimilarity',pch=19)
abline(lm(dist(as.matrix(ws$wses))~vegdist(cbind(wco,ds=rep(0.0001,nrow(wco))))))
                                        #ch~ses
plot(dist(as.matrix(ws$wses))~vegdist(wch),xlab='Checker Dissimilarity',
     ylab='SES Dissimilarity',pch=19)
abline(lm(dist(as.matrix(ws$wses))~vegdist(wch)))
## Araujo network
wan <- coNet(do.call(rbind,wq))
wan <- log(wan+1)
wan.col <- coNet(do.call(rbind,wq),showsign=TRUE)
wan.vs <- (apply(wc,2,sum))
wan.vs <- wan.vs/max(wan.vs)+1.5
wan.vc <- c(wom.[[2]],3)
wan.vc <- wan.vc[match(names(wan.vs),names(wan.vc))]
wan.vc <- cbp[wan.vc]
wan.vc[is.na(wan.vc)] <- 'grey'
names(wan.vc)[wan.vc=='grey'] <- 'Pu'
coord <- coord[match(rownames(wan),rownames(coord)),]
gplot(wan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='lightgrey',
      vertex.cex=wan.vs,edge.lwd=wan)
text(coord,rownames(wan))

###average network
waan <- lapply(wq,coNet)
waan <- Reduce('+',waan)
waan <- log(waan+1)
coord <- coord[match(rownames(waan),rownames(coord)),]
gplot(waan,displaylabels=FALSE,gmode='graph',coord=coord,vertex.col='lightgrey',
      vertex.cex=wan.vs,edge.lwd=waan)
text(coord,rownames(waan))
sum(abs(sign(wan)-sign(waan)))


###Eigen Centrality
library(igraph)
wnec <- centralization.evcent(graph.adjacency(wan,weighted=TRUE))$vector
plot(graph.adjacency(wan,weighted=TRUE),vertex.size=wnec*10)
names(wnec) <- rownames(wan)
round(wnec,2)
plot(wnec[match(names(onec),names(wnec))]~onec,pch=19,xlab='Eigen Centrality (ONC)',ylab='Eigen Centrality (Uintah)',xlim=c(0,1))
abline(lm(wnec[match(names(onec),names(wnec))]~onec),lty=1)
cor.test(wnec[match(names(onec),names(wnec))],onec,method='p')

###Abundance ~ centrality
cor.test(wnec,log(wan.vs+1),method='p')
cor.test(onec,log(oan.vs+1),method='p')
                                        #
plot(onec~I(log(oan.vs+1)/max(log(oan.vs+1))),pch=19,xlab='Relative log(Abundance)',ylab='Eigen Centrality')
points(wnec~I(log(wan.vs+1)/max(log(wan.vs+1))),pch=19,col='darkgrey')
abline(lm(onec~I(log(oan.vs+1)/max(log(oan.vs+1)))))
abline(lm(wnec~I(log(wan.vs+1)/max(log(wan.vs+1)))),lty=2)
legend('bottomright',legend=c('Garden','','Natural',''),pch=c(19,31,19,31),col=c(1,1,'darkgrey',1),lty=c(0,1,0,2))
