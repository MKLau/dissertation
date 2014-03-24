###LCN: Wild Stand Analyses
###MKLau
###21Mar2014

###Meta
###Site = Uintah, UT
###Study area = 225 * 463 = 104,175 m2 = 0.104175 km2

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

wmods <- oecosimu(apply(wc[,colnames(wc)!='Pu'],2,function(x) x/max(x)),mm,method='r1',burnin=100,thin=10,nsimul=1000)
save(wmods,'../data/lcn_mods_wild.rda')
