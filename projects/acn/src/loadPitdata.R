### script for loading data for the
### ACN study
### Pit data collected 2012

source('../src/helpers.R')
pit <- read.csv('../data/arth_cooc_PIT_Lau.csv')
pit[is.na(pit)] <- 0

### separating out "fungal" necrosis
necrosis <- pit$fungal
pit <- pit[,colnames(pit) != 'fungal']

### combine pb 
pb <- pit$pb.upper + pit$pb.lower + pit$pb.woody 
pb.pred <- pit$pb.pred + pit$pb.woody.pred + pit$pb.hole
pb.abort <- pit$pb.abort
pit <- pit[,!(grepl('pb',colnames(pit)))]
pit <- data.frame(pit,pb.abort,pb.pred,pb)
pit <- data.frame(pit[,1:6],pit[,ncol(pit):7])

### separating into trees
tree.info <- pit[,c(1,2,3,6)]
tree.arth <- pit[,7:ncol(pit)]
tree.arth <- split(tree.arth,paste(tree.info[,1],
                                   tree.info[,2],
                                   tree.info[,3]))
tree.info <- tree.info[!(duplicated(tree.info)),]

### coerce into numeric matrices
arth.mats <- lapply(tree.arth,data.frame)
arth.mats <- lapply(arth.mats,as.matrix)

### species totals for each tree
spp.tot <- do.call(rbind,lapply(tree.arth,function(x) apply(x,2,sum)))
