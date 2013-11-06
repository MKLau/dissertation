###Arthropod Co-occurrence Networks
###Null modeling to get ses and p-values

########################################################################
print('Loading packages')
source('../../lichen_coo/src/seenetR.R')
library(vegan)
library(pbapply)
########################################################################
print('Loading data..')
pit <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc_PIT_Lau.csv')
pit[is.na(pit)] <- 0
pit.com <- pit[,-c(1,2,3,4,5)]
pit.l <- split(pit.com,paste(pit$tree,pit$geno))
########################################################################
print('Tree level co-occurrence')
obs.cs <- unlist(lapply(pit.l,cscore))
acn.sim <- pblapply(pit.l,function(x) if (sum(sign(apply(x,2,sum)))>1){nullCom(x)}else{NA})
acn.cs <- pblapply(acn.sim,function(x) if (any(is.na(x[[1]]))){NA}else{lapply(x,cscore)})
acn.cs <- pblapply(acn.cs,unlist)
acn.ses <- obs.cs*0
acn.p <- obs.cs*0
for (i in 1:length(acn.ses)){
  print(i)
  acn.ses[i] <- (obs.cs[i] - mean(acn.cs[[i]])) / sd(acn.cs[[i]])
  acn.p[i] <- length(acn.cs[[i]][acn.cs[[i]]<=obs.cs[i]])/length(acn.cs[[i]])
}
print('Writing output')
dput(acn.cs,'../data/acn_tree_cs_nb.Rdata')
dput(acn.ses,'../data/acn_tree_ses_nb.Rdata')
dput(acn.p,'../data/acn_tree_pval_nb.Rdata')
