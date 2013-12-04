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
pit$tree <- sub('\\.0','\\.',pit$tree)
pit <- pit[pit$geno!='1007',]
pit[is.na(pit)] <- 0
                                        #merge categories
                                        #pemphigus mergers
pit$pb.upper <- pit$pb.upper + pit$pb.woody
pb <- pit$pb.upper + pit$pb.lower
pit <- cbind(pit,pb=pb)
pit$pb.pred <- pit$pb.pred + pit$pb.hole + pit$pb.woody.pred
pit <- pit[,colnames(pit)!='pb.woody'&colnames(pit)!='pb.woody.pred'&colnames(pit)!='pb.hole'&colnames(pit)!='mite'&colnames(pit)!='pb.upper'&colnames(pit)!='pb.lower']
                                        #remove species with less than 17 observations
pit.com <- pit[,-1:-6]
pit.com <- pit.com[,apply(pit.com,2,sum)>17]
                                        #separate live and senescing leaves
liv <- pit.com[pit[,1]=='live',]
pit.l <- split(liv,paste(pit$tree[pit[,1]=='live'],pit$geno[pit[,1]=='live']))
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
dput(acn.cs,'../data/acn_cs_live.Rdata')
dput(acn.ses,'../data/acn_ses_live.Rdata')
dput(acn.p,'../data/acn_pval_live.Rdata')
