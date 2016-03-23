### script for loading data for the
### ACN study
### Pit data collected 2012

source('helpers.R')
pit <- read.csv('../data/arth_cooc_PIT_Lau.csv')
pit[is.na(pit)] <- 0

### separating out "fungal" necrosis
necrosis <- pit$fungal
pit <- pit[,colnames(pit) != 'fungal']

### separating into trees
tree.info <- pit[,1:6]
tree.arth <- pit[,7:ncol(pit)]
tree.arth <- split(tree.arth,paste(tree.info[,1],tree.info[,2],tree.info[,3]))

### generate network models
arth.sub <- lapply(tree.arth,zeroCol,n=10)
do.call(rbind,lapply(arth.sub,function(x) sum(sign(apply(x,2,sum)))))
arth.cor <- lapply(arth.sub,cor,method='pearson')
