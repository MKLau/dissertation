###Genotype sensitivity

binned.species <- TRUE

library(sna)
source('../../art_coo/src/helper_func.R')
source('../../lichen_coo/src/seenetR.R')

###Prelim data from ONC
## onc <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc.csv')
## onc.key <- read.csv('~/projects/dissertation/projects/art_coo/data/key.csv')
###Data from Pit

pit <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc_PIT_Lau.csv')
                                        #genotype remove 1007
pit <- pit[pit$geno!='1007',]

###Data checks
                                        #na to zeros
pit[is.na(pit)] <- 0
                                        #merge categories
pit.com <- pit[,-c(1,2,3,4,5)]
if (binned.species){
  pb <- apply(pit.com[,1:7],1,sum)
  pit.com <- cbind(pb,pit.com[,-1:-7])
  pit.com <- cbind(pit.com[,-c(2,3)],chew=apply(pit.com[,2:3],1,sum))
  pit.com <- pit.com[,apply(pit.com,2,sum)>1]
  pit.com$pb[pit.com$pb!=0] <- 1
  pit.com$chew[pit.com$chew!=0] <- 1
  pit.com <- as.matrix(pit.com)
}else{}
                                        #
pit.l <- split(pit.com,paste(pit$tree,pit$geno))
geno <- as.character(unlist(sapply(names(pit.l),function(x) strsplit(x,split=' ')[[1]][2])))
                                        #
pit.dnet <- dep.net(pit.com)
                                        #Genotype sensitivity
                                        #Random tree removal
net.all <- pit.dnet
ng <- 5
for (k in 1:ng){
net.rnd <- list()
n <- round(mean(table(geno)),0)*k
perm <- 5000
for (i in 1:perm){
                                        #randomly remove trees across genotypes
  rnd.tree <- sample(pit$tree,n)
  geno.pc <- table(pit$gen[match(rnd.tree,pit$tree)])/sum(table(pit$gen[match(rnd.tree,pit$tree)]))
  while (any(geno.pc >= 0.5)){
      rnd.tree <- sample(pit$tree,n)
      geno.pc <- table(pit$gen[match(rnd.tree,pit$tree)])/sum(table(pit$gen[match(rnd.tree,pit$tree)]))
  }
  net.rnd[[i]] <- dep.net(pit.com[-match(rnd.tree,pit$tree),])
}
dput(net.rnd,file=paste('../data/acn_net_rnd_',k,'.Rdata',sep=''))
}
