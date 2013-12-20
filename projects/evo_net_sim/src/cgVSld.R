###Comparing the cgSim output to the lonsdorf octave output
source('~/projects/dissertation/projects/lichen_coo/src/seenetR.R')
library(vegan)
library(rgl)
trees <- read.csv('../data/trees.txt')
geno <- letters[trees[,1]]
                                        #lonsdorf output
ld <- list()
for (i in 1:length(paste('../data/lonsdorf_out/',dir('../data/lonsdorf_out/'),sep=''))){
  ld[[i]] <- as.matrix(read.table(paste('../data/lonsdorf_out/',dir('../data/lonsdorf_out/'),sep='')[i]))
}
ld.names <- dir('../data/lonsdorf_out/')
ld.num <- as.numeric(sapply(ld.names,function(x) strsplit(x,split='_')[[1]][[2]]))
ld <- ld[order(ld.num)]
                                        #comsim output
                                        #reps*GG*(y-1),reps*(z-1),RR
cg <- dget(file='../data/cg_out/cg_output')
cg.l <- list()
l <- 1
for (i in 1:length(cg)){
  for (j in 1:length(cg[[i]])){
    for (k in 1:length(cg[[i]][[j]])){
      print(l)
      cg.l[[l]] <- cg[[i]][[j]][[k]]
      names(cg.l)[l] <- names(cg[[i]][[j]])[k]
      l <- l + 1
    }
  }
}
cg.perm <- lapply(cg.l,function(x,f) adonis(x~f,permutations = 3)$aov.tab,f=geno)
cg.r2 <- unlist(lapply(cg.perm,function(x) x[1,5]))
cg.yy <- as.numeric(unlist(sapply(names(cg.l),function(x) strsplit(x,split='_')[[1]][2])))
cg.gg <- as.numeric(unlist(sapply(names(cg.l),function(x) strsplit(x,split='_')[[1]][3])))
pairs(cbind(cg.yy,cg.gg,cg.r2))
plot3d(cbind(y=cg.yy,z=cg.gg,cg.r2),col=terrain.colors(7)[round(cg.r2*7,0)])
plot3d(cbind(y=(((cg.yy-1)*VeN-VeN*(cg.yy-1)/2)),z=(cg.gg-1)*(15),cg.r2),col=terrain.colors(7)[round(cg.r2*7,0)])

###Bipartite networks
library(bipartite)
library(vegan)
bpn.resolve <- function(x,alpha=0.05){
  p <- t.test(x,alternative='great')$p.value
  if (p<alpha){return(mean(x))}else{return(0)}
}
test <- cg[[10]][[5]][[8]]
bpn <- array(NA,dim=c(length(unique(geno)),ncol(test)))
for (i in 1:length(unique(geno))){
  for (j in 1:ncol(bpn)){
    bpn[i,j] <- bpn.resolve(test[geno==unique(geno)[i],j],alpha=0.05/length(unique(geno)))
  }
}
plotweb(bpn)
visweb(bpn)
nestedchecker(bpn)
nest.r0 <- oecosimu(bpn,nestedtemp,method='r0',nsimul=1000,burnin=100)
nest.r0
cg.nest.10.5.7 <- nest.r0

###Unipartite networks
test <- cg.l
thresh <- 5
for (i in 1:length(test)){
  test[[i]][test[[i]]<=thresh] <- 0
}
nets <- lapply(test,CoNetwork,plot.net=FALSE)
par(mfcol=c(2,3))
for (i in (length(nets)-2):length(nets)){
  e.col <- nets[[i]]
  e.col[nets[[i]]!=0&nets[[i]]<0.5] <- 'black'
  e.col[nets[[i]]!=0&nets[[i]]>0.5] <- 'red'
  gplot(abs(nets[[i]]),displaylabels=TRUE,gmode='graph',edge.col=e.col)
  hist(nets[[i]],main='')
}

par(mfcol=c(2,3))
net.choose <- (1:length(cg.r2))[cg.r2>0.5]
net.choose <- sample(net.choose,3)
for (i in net.choose){
  e.col <- nets[[i]]
  e.col[nets[[i]]!=0&nets[[i]]<0.5] <- 'black'
  e.col[nets[[i]]!=0&nets[[i]]>0.5] <- 'red'
  gplot(abs(nets[[i]]),displaylabels=TRUE,gmode='graph',edge.col=e.col)
  hist(nets[[i]],main=round(cg.r2[i],2))
}
net.deg <- unlist(lapply(nets,function(x) sum(abs(sign(x)))/2))
pairs(cbind('Network Size'=net.deg,'selection intensity'=cg.gg,'environmental variation'=cg.yy,'R2 PerMANOVA'=cg.r2))
