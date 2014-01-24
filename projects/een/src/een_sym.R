###EEN = evolution on ecological networks
###Does genotypic diversity produce robustness through nestedness?
###Simulation experiments
###Initiated 23 Jan 2014

### Simulate networks with varying levels of genotypic effects 
## Use percent genetic variance equivalent to published community heritability estimates 
## Whitham et al. 2013 (H2C ~ 65% for arthropods)
## Tree clonal phenotype heritability (52% from Lamit's Roughness, similar to tannins from NRG2006)

output.loc <- '../results/een_exp_sym/'
unlink(output.loc,recursive=TRUE)
dir.create(output.loc)
round.sim <- TRUE

library(ComGenR)
tree.gpm <- gpmTrees(reps=10)
#make H2B=55%
trees <- simTrees(tree.gpm,VeT=7.5);h2b <- getH2C(trees,tree.gpm[,1],method='nmds')[2]
while (h2b<0.55|h2b>0.60){trees <- simTrees(tree.gpm,VeT=7.5);h2b <- getH2C(trees,tree.gpm[,1],method='nmds')[2]}
insects <- simSpp()
                                        #write trees and insects out
write.csv(tree.gpm,file='../data/tree_gpm.csv',row.names=FALSE)
write.csv(trees,file='../data/trees.csv',row.names=FALSE)
write.csv(insects,file='../data/insects.csv',row.names=FALSE)
                                        #generate communities
cg.sim <- out <- list()
for (i in 1:9){
  for (j in 1:10){
    out[[j]] <- cgSim(trees=trees,insects=insects,z=i,Ve=0.1,VeN=3,k.asym=FALSE)
    print(paste(i,j))
  }
  cg.sim[[i]] <- out
}
threshold <- 5
com.sim <- list()
for (i in 1:length(cg.sim)){
  for (j in 1:length(cg.sim[[1]])){
    out[[j]] <- cg.sim[[i]][[j]]$art_pop
    out[[j]][out[[j]]<threshold] <- 0
    if (round.sim){out[[j]] <- round(out[[j]],0)}else{}
  }
  com.sim[[i]] <- out
}

###For an appendix, develop the use of permanova R2, inlcuding:
##Calculation using Anderson method
##Comparison to NMDS method in Shuster 2006
g <- tree.gpm[,1]
h2c.shuster <- lapply(com.sim,function(x) lapply(x,function(x,g) getH2C(x,g,method='nmds'), g=g))
h2c.perm <- pblapply(com.sim,function(x) lapply(x,function(x,g) getH2C(x,g,method='permanova'), g=g))
write.csv(h2c.shuster,file=paste(output.loc,'h2c_shuster.csv',sep=''),row.names=FALSE)
write.csv(h2c.perm,file=paste(output.loc,'h2c_perm.csv',sep=''),row.names=FALSE)

###Measure nestedness and examine correlation between genetic variance and nestedness

x <- com.sim
cargo <- list()
out <- list()
for (i in 1:length(x)){
  for (j in 1:length(x[[i]])){
    if (any(x[[i]][[j]]==0)){cargo[[j]] <- nestedtemp(x[[i]][[j]])$statistic}else{cargo[[j]] <- 0}
    print(paste(i,j))
  }
  out[[i]] <- cargo
}

nest.sim <- cbind(z=as.numeric(gl(length(out),length(out[[1]]))),nestedtemp=as.numeric(unlist(out)))
write.csv(nest.sim,file=paste(output.loc,'nest_sim.csv',sep=''),row.names=FALSE)

###Conduct Removal experiments (random vs targeted removal):

###Threshold of species loss = number of trees removed until a species goes locally extinct
rm.nits <- 5000

cargo <- list()
##Random removal
rnd.rm <- list()
for (i in 1:length(com.sim)){
  for (j in 1:length(com.sim[[i]])){
    cargo[[j]] <- rmTrees(com.sim[[i]][[j]],method='random',nits=rm.nits)
    names(cargo)[j] <- paste(i,j)
  }
  rnd.rm[[i]] <- unlist(cargo)
  print(paste(i,j))
}
rnd.rm <- unlist(rnd.rm)
##Degree removal
deg.rm <- list()
for (i in 1:length(com.sim)){
  for (j in 1:length(com.sim[[i]])){
    cargo[[j]] <- rmTrees(com.sim[[i]][[j]],method='degree',nits=rm.nits)
    names(cargo)[j] <- paste(i,j)
  }
  deg.rm[[i]] <- unlist(cargo)
  print(paste(i,j))
}
deg.rm <- unlist(deg.rm)
##Type removal
typ.rm <- list()
for (i in 1:length(com.sim)){
  for (j in 1:length(com.sim[[i]])){
    cargo[[j]] <- rmTrees(com.sim[[i]][[j]],method='type',type=g,nits=rm.nits)
    names(cargo)[j] <- paste(i,j)
  }
  typ.rm[[i]] <- unlist(cargo)
  print(paste(i,j))
}
typ.rm <- unlist(typ.rm)
write.csv(cbind(rnd.rm,deg.rm,typ.rm),file=paste(output.loc,'removal_results.csv',sep=''))
