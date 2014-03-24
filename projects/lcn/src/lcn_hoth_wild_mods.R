###LCN: Wild Stand Analyses
###MKLau
###21Mar2014

###Meta
###Site = Uintah, UT
###Study area = 225 * 463 = 104,175 m2 = 0.104175 km2
library(bipartite)
library(methods)
mm <- function(x){slot(computeModules(x),'likelihood')}
                                        #
x <- read.csv('~/projects/dissertation/projects/lcn/data/lco_Apr2012.csv')
                                        #remove notes
x <- x[,colnames(x)!='NOTES.']
x <- x[,colnames(x)!='dead']
                                        #
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
x$tree <- factor(as.character(x$tree))
                                        #condense species
                                        #lecanora, there can be only one!
lec.sp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
                                        #no physcioids!
                                        #phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x,lec=lec.sp)
x <- x[,-c(6,8,10,18)]
x <- x[,colnames(x)!='physcioid']
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
colnames(x)[5:ncol(x)] <- c('Xg','Cs', 'Xm', 'fgb', 'Rs', 'Pm' ,'Pa', 'Pu','Ch','Ls')
x <- x[colnames(x)!='fgb']
x.q <- split(x,quads)
wild.com <- split(x,x$tree)
wild.com <- do.call(rbind,lapply(wild.com,function(x) apply(x[,-1:-4],2,sum)))

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
wc <- wild.com
wmods <- oecosimu(apply(wc[,colnames(wc)!='Pu'],2,function(x) x/max(x)),mm,method='r1',burnin=100,thin=10,nsimul=1)
save(wmods,'../data/lcn_mods_wild.rda')
