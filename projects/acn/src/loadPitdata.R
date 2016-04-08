### script for loading data for the
### ACN study
### Pit data collected 2012

source('../src/helpers.R')
library(sna)
library(vegan)
library(ecodist)

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

### tree level networks for arthropods
tree.nets <- lapply(arth.mats,coNets)


## Genotype average network 
## within each 
liv.nets <- tree.nets[tree.info[,1] == 'live']
sen.nets <- tree.nets[tree.info[,1] == 'sen']
for (i in 1:length(liv.nets)){
    diag(liv.nets[[i]]) <- 0
    diag(sen.nets[[i]]) <- 0
}

liv.evc <- unlist(lapply(liv.nets,function(x) centralization(x,FUN='evcent')))
sen.evc <- unlist(lapply(sen.nets,function(x) centralization(x,FUN='evcent'))
)

liv.dc <- unlist(lapply(liv.nets,function(x) centralization(x,FUN='degree')))
sen.dc <- unlist(lapply(sen.nets,function(x) centralization(x,FUN='degree')))

### species eigen centralities
liv.spc <- list()
for (i in 1:length(liv.nets)){
        liv.spc[[i]] <- evcent(liv.nets[[i]],gmode='graph')
}
liv.spc <- do.call(rbind,liv.spc)
liv.spc[is.na(liv.spc)] <- 0
colnames(liv.spc) <- colnames(liv.nets[[1]])

sen.spc <- list()
for (i in 1:length(sen.nets)){
        sen.spc[[i]] <- evcent(sen.nets[[i]],gmode='graph')
}
sen.spc <- do.call(rbind,sen.spc)
sen.spc[is.na(sen.spc)] <- 0
colnames(sen.spc) <- colnames(sen.nets[[1]])

### edge degree

liv.ord <- unlist(lapply(liv.nets,ord))
sen.ord <- unlist(lapply(sen.nets,ord))

