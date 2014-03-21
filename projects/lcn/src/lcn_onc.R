###LCN: ONC Garden Analyses
###MKLau
###21Mar2014

source('~/projects/dissertation/projects/lcn/src/lcn_load_onc.R')
source('~/projects/packages/cooc/src/cooc.R')

###Data objects:
###oc = co-occurrences summed across all cells for each tree
###oq = co-occurrence matrices separated out for each tree
###os = matrix of ses and related stats for trees
###wco = co-occurrence counts
###wch = checker counts
###prb = percent rough bark (averaged between the upper and lower)
###Data notes:
###No physciods
###Lecanoras merged

## Height composition
### Previous tests show that this is not significant

## Overall Co-occurrence Patterns
#wco.all <- oecosimu(do.call(rbind,wq),cs,method='r1',burnin=100,thin=10,nsimul=5000)
oecosimu(oc,cs,method='r1',burnin=100,thin=10,nsimul=100)

## oecosimu(wc,cs,method='r1',burnin=100,thin=10,nsimul=5000)
##           statistic       z   mean   2.5%    50%  97.5%  Pr(sim.)    
## statistic    1.3556 -2.9207 4.5655 2.3778 4.5556 6.7339 0.0005999

## Nestedness, Modularity, Centrality
### oecosimu(wc,nested,method='r1',burnin=100,thin=10,nsimul=5000)
##             statistic       z   mean   2.5%    50%  97.5% Pr(sim.)   
## binmatnest2    8.6663 -2.5715 18.592 11.589 18.428 26.445   0.0018 **
## Only did 10
## oecosimu(wc,mm,method='r1',burnin=100,thin=10,nsimul=10)
##           statistic       z    mean    2.5%     50%  97.5% Pr(sim.)
## statistic   0.18652 -1.6232 0.20923 0.18708 0.21407 0.2263   0.2727
## mods <- computeModules(wc)

plotweb(wc[order(apply(wc,1,sum),decreasing=TRUE),order(apply(wc,2,sum),decreasing=TRUE)],text.rot=90,method='normal')

## SES ~ Roughness
## prb         -0.03937    0.01809  -2.177   0.0502 .
summary(lm(ws$wses~prb))

shapiro.test(residuals(lm(ws$wses~prb)))
plot(ws$wses~prb,xlab='Percent Rough Bark',ylab='SES',pch=19)
abline((lm(ws$wses~prb)))

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
wan <- newAraujo(do.call(rbind,wq)[,c(-4,-8)])
wan.col <- newAraujo(do.call(rbind,wq)[,c(-4,-8)],showsign=TRUE)
gplot(,displaylabels=TRUE,gmode='graph')
