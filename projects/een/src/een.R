###EEN = evolution on ecological networks
###Does genotypic diversity produce robustness through nestedness?
###Simulation experiments
###Initiated 23 Jan 2014

### Simulate networks with varying levels of genotypic effects 
## Use percent genetic variance equivalent to published community heritability estimates 
## Whitham et al. 2013 (H2C ~ 65% for arthropods)
## Tree clonal phenotype heritability (52% from Lamit's Roughness, similar to tannins from NRG2006)

library(ComGenR)

tree.gpm <- gpmTrees(reps=10)
#make H2B=55%
trees <- simTrees(tree.gpm,VeT=7.5);h2b <- getH2C(trees,tree.gpm[,1])[2]
while (h2b<0.55|h2b>0.60){trees <- simTrees(tree.gpm,VeT=7.5);h2b <- getH2C(trees,tree.gpm[,1])[2]}
insects <- simSpp()
max.k <- 100
min.k <- 5
p.k <- 1-ppois(1:nrow(insects),5)
k.vector <- (p.k*max.k) + ((1-p.k)*min.k)
k.vector <- round(k.vector,0)
                                        #generate communities
cg.sim <- out <- list()
for (i in 1:9){
  for (j in 1:10){
    out[[j]] <- cgSim(trees=trees,insects=insects,z=i,Ve=0.1,VeN=3,K=k.vector)
    print(paste(i,j))
  }
  cg.sim[[i]] <- out
}
threshold <- min.k
com.sim <- list()
for (i in 1:length(cg.sim)){
  for (j in 1:length(cg.sim[[1]])){
    out[[j]] <- cg.sim[[i]][[j]]$art_pop
    out[[j]][out[[j]]<threshold] <- 0
  }
  com.sim[[i]] <- out
}


###For an appendix, develop the use of permanova R2, inlcuding:
##Calculation using Anderson method
##Comparison to NMDS method in Shuster 2006

g <- tree.gpm[,1]
##h2c.shuster <- lapply(com.sim,function(x) lapply(x,function(x,g) getH2C(x,g,method='nmds'), g=g))
##h2c.perm <- pblapply(com.sim,function(x) lapply(x,function(x,g) getH2C(x,g,method='permanova'), g=g))
h2c.shuster <- read.csv(file='../results/h2c_shuster.csv')
h2c.perm <- read.csv(file='../results/h2c_perm.csv')

plot(h2c.shuster$H2C~h2c.shuster$z,pch=19,col='grey',ylim=c(-0.1,0.65),xlab='Selection Intensity',ylab='H2C')
points(h2c.perm$H2C~h2c.perm$z,pch=19,col='black',ylim=c(0,1))

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

## nest.sim <- cbind(nestedtemp=as.numeric(unlist(out)),z=as.numeric(gl(length(out),length(out[[1]]))))
## write.csv(nest.sim,file='../results/nest_sim.csv',row.names=FALSE)
nest.sim <- read.csv(file='../results/nest_sim.csv')
plot(nest.sim[,1]~nest.sim[,2],pch=1,xlab='Selection Intensity',ylab='Nested Temperature',cex=0.75)

###Conduct Removal experiments (random vs targeted removal):

###Threshold of species loss = number of trees removed until a species goes locally extinct

##Random removal
rnd.rm <- list()
for (i in 1:length(com.sim)){
  j <- sample(1:length(com.sim[[i]]),1)
  rnd.rm[[i]] <- rmTrees(com.sim[[i]][[j]])
  print(paste(i,j))
  names(rnd.rm)[i] <- paste(i,j)
}
rnd.rm <- unlist(rnd.rm)
avg.temp <- tapply(nest.sim[,1],nest.sim[,2],mean)
par(mfrow=c(1,2))
plot(rnd.rm~I(1:9),xlab='Selection Intensity',ylab='Percent Removed')
plot(rnd.rm~avg.temp,xlab='Avg. Nestedness',ylab='Percent Removed')

step

##Random removal
##Genotype removal based on phenotpic similarity
##Centralized species removal
