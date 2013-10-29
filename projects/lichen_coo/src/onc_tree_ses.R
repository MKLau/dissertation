##Generate the stand level ses values
source('./seenetR.R')
library(vegan)
library(pbapply)
print('Loading data')
garden.data <- read.csv('../data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
					#remove RL6 and N1.31
print('remove genotype RL6 and tree N1.31')
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Geno!='N1.31',]
                                        #separate onc
print('ONC analyses')
onc <- garden.data[garden=='onc',]
###Merge species groups
print('Merge physcia species')
Phy <- apply(onc[,12:14],1,sum) #make phy out of pmel,pads,pund
onc. <- onc[,-12:-14]
onc <- data.frame(onc.,Phy)
					#tree level
print('Tree level co-occurrence')
onc.q <- split(onc,paste(onc[,1],onc[,3],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
obs.cs <- unlist(lapply(onc.q,cscore))
onc.sim <- pblapply(onc.q,function(x) if (sum(sign(apply(x,2,sum)))>1){nullCom(x)}else{NA})
onc.cs <- pblapply(onc.sim,function(x) if (any(is.na(x[[1]]))){NA}else{lapply(x,cscore)})
onc.cs <- pblapply(onc.cs,unlist)
onc.ses <- obs.cs*0
onc.p <- obs.cs*0
for (i in 1:length(onc.ses)){
  print(i)
  onc.ses[i] <- (obs.cs[i] - mean(onc.cs[[i]])) / sd(onc.cs[[i]])
  onc.p[i] <- length(onc.cs[[i]][onc.cs[[i]]<=obs.cs[i]])/length(onc.cs[[i]])
}
print('Writing to: ../data/onc_tree_ses.Rdata')
dput(onc.cs,'../data/onc_tree_cs.Rdata')
dput(onc.ses,'../data/onc_tree_ses.Rdata')
dput(onc.p,'../data/onc_tree_pval.Rdata')

