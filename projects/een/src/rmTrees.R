###rmTrees
###Tree removal simulator
###Three options for removal: 1) random, 2) degree, 3) (geno-, pheno-)type

rmTrees <- function(x='network',method=c('random','degree','type'),type='grouping'){
  if (method[1]=='random'){
    rm.x <- x[,apply(x,2,function(x) sum(x))!=0]
    rm.x[rm.x!=0] <- 1
    itr <- 0
    spp.deg <- apply(rm.x,2,sum)
    live.trees <- sign(apply(rm.x,1,sum))
    rm.prob <- live.trees/length(live.trees[live.trees!=0])
    while(all(spp.deg!=0)){
      live.trees <- sign(apply(rm.x,1,sum))
      rm.prob <- live.trees/length(live.trees[live.trees!=0])
      rm.tree <- sample((1:nrow(rm.x)),1,prob=rm.prob)
      rm.x[rm.tree,] <- 0
      spp.deg <- apply(rm.x,2,sum)
      itr <- itr + 1
      image(t(rm.x))
    }
    out <- (itr/nrow(x))
  }else if (method[1]=='degree'){
    rm.x <- x[,apply(x,2,function(x) sum(x))!=0]
    rm.x[rm.x!=0] <- 1
    itr <- 0
    spp.deg <- apply(rm.x,2,sum)
    deg.trees <- apply(rm.x,1,sum)
    rm.prob <- deg.trees/max(deg.trees)
    while(all(spp.deg!=0)){
      deg.trees <- apply(rm.x,1,sum)
      rm.prob <- deg.trees/max(deg.trees)
      rm.tree <- sample((1:nrow(rm.x)),1,prob=rm.prob)
      rm.x[rm.tree,] <- 0
      spp.deg <- apply(rm.x,2,sum)
      itr <- itr + 1
      image(t(rm.x))
    }
    out <- (itr/nrow(x))
  }else if (method[1]=='type'){
    if (type=='grouping'){warning('Please provide a grouping vector.')}
    rm.x <- x[,apply(x,2,function(x) sum(x))!=0]
    rm.x[rm.x!=0] <- 1
    itr <- 0
    spp.deg <- apply(rm.x,2,sum)
    rnd.type <- sample(unique(type),1)
    rnd.tree <- (1:length(type))[type==rnd.type][1]
    rm.prob <- as.matrix(dist(type))[rnd.tree,]
    rm.prob <- rm.prob/max(rm.prob)
    while(all(spp.deg!=0)){
      live.trees <- apply(rm.x,1,function(x) sign(sum(x)))
      rnd.type <- sample(unique(type),1)
      rnd.tree <- (1:length(type))[type==rnd.type][1]
      rm.prob <- as.matrix(dist(type))[rnd.tree,]
      rm.prob <- rm.prob/max(rm.prob)
      rm.prob[live.trees==0] <- 0
      rm.tree <- sample((1:nrow(rm.x)),1,prob=rm.prob)
      rm.x[rm.tree,] <- 0
      spp.deg <- apply(rm.x,2,sum)
      itr <- itr + 1
      image(t(rm.x))
    }
    out <- (itr/nrow(x))
  }else{warning('Unknown simulation method.')}
  return(out)
}
