###Functions for the SEEnetR package
###Spatial, Ecological and Evolutionary Networks in R
###25 Sep 2013

require(vegan)
require(igraph)

###Spatial

null.prune <- function(a,b,std=TRUE){

###Method for producing network models from co-occurrence data.
###Developed in Araujo et al. 2011 Ecography Using species co-occurrence networks to assess the impacts of climate change
###Coded 18 Sep 2013 MKLau

## E is an environmental lattice with rows = N and cols = M and has size A = N*M.
## A = total number of sites
## Na = number of sites with species a
                                        #determine number of sites
A <- length(a)
                                        #NULL probabilities
## Assuming no interactions among species
## Probability of observing species a
## Pa <- Na/A
## Probability of observing both species a and b simultaneously
## Pa&b = Pa * Pb
## Probability of observing either species a and b
## Pa+b = Pa(1-Pb) + (1-Pa)Pb = Pa + Pb - 2PaPb 
## Probability of observing neither species a nor b
## P!a&!b = (1-Pb)(1-Pa) = 1 - Pa - Pb + PaPb
## and it follows that 
## Pa&b + Pa+b + P!a&!b = 1

pa <- sum(a)/A
pb <- sum(b)/A
paANDb <- pa * pb
paORb <- pa + pb - 2*pa*pb
pNOTaORb <- 1- pa - pb + pa*pb

## Three states, i = ab, ii a+b = (a&!b + !a&b), iii = !a!b
                                        #expected frequencies
## E(i) = APa&b
## E(ii) = APa+b
## E(iii) = AP!a&!b

Ei <- A*paANDb
Eii <- A*paORb
Eiii <- A*pNOTaORb
                                        #variance of frequencies
## Var(i) = APa&b(1-Pa&b)
## Var(ii) = APa+b(1-APa+b)
## Var(iii) = AP!a&!b(1-AP!a&!b)

Vi <- A*paANDb*(1-paANDb)
Vii <- A*paORb*(1-paORb)
Viii <- A*pNOTaORb*(1-pNOTaORb)
                                        
## Critical interval
## Critical threshold (+ = attraction and - = repulstion)
## CI.u = E(i) + 2*Var(i)^2 = attractive overlap
## CI.l = E(i) - 2*Var(i)^2 = repulsive overlap

ci.u <- Ei + 2*Vi^(1/2)
ci.l <- Ei - 2*Vi^(1/2)
                                        #Oi = observed a AND b
Oi <- length(a[(a+b)==2])
                                        #Determine output (Oi = significant, 0 = not)
if (Oi > ci.u | Oi < ci.l){
  if (std){
    out <- (Oi - Ei) / sqrt(Vi)
  }else{
    out <- Oi
  }

}else{
  out <- 0
}
return(out)
}

co.net <- function(x='species in cols',diag.zero=TRUE,std=TRUE){
  x[x!=0] <- 1
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  rownames(out) <- colnames(out) <- colnames(x)
  for (j in 1:ncol(x)){
    for (k in 1:ncol(x)){
      out[j,k] <- null.prune(x[,j],x[,k],std=std)
    }
  }
  if (diag.zero){diag(out) <- 0}else{}
  return(out)
}


##Species Dependence
calcDepend <- function(a,b){
  return(length(a[(a+b)==2])/sum(a)) #intersection of a with b divided by the total of a
}

###Dependency Network
dep.net <- function(x='species in cols',zero.na=TRUE,prune=TRUE,diag.zero=TRUE,pos=TRUE){
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (pos){
        out[i,j] <- calcDepend(x[,i],x[,j])
      }else{
        out[i,j] <- negDepend(x[,i],x[,j])
      }
    }
  }
  if (prune){
    out.rm <- co.net(x,diag.zero=diag.zero)
    out[out.rm==0] <- 0
  }else{}
  if (diag.zero){diag(out) <- 0}
  rownames(out) <- colnames(out) <- colnames(x)
  if (zero.na){out[is.na(out)] <- 0}
  return(out)
}

percThreshold <- function(x='network matrix',step.size=0.01){
  no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
  step <- 1
  while (no.c==1){
    x[x<=(step*step.size)] <- 0
    no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
    step <- step + 1
  }
  out <- list(threshold=((step-1)*step.size),isolated.nodes=colnames(x)[clusters(graph.adjacency(x,weighted=TRUE))$membership==2])
  return(out)
}


CoNetwork <- function(x,plot.net=TRUE,scalar=3,min.vsize=0.1){
###Runs all steps of the process for modeling
###Co-occurrence networks described by Araujo et al. 2011.
###It depends on the seenetR.R script which contains both the
###Araujo functions and related co-occurrence null modeling
###functions.

###Inputs: 
#x = matrix of co-occurrence with species in columns
#plot.net = logical. Should the network be plotted?
#scalar = scales the size of all vertices
#min.vsize = sets the minimum size for vertices

#Step 1. Calculate a Bray-Curtis distance matrix

bc.d <- as.matrix(vegdist(t(x)))

#Step 2. Prune distance matrix based on co-occurrence probabilities

prune <- co.net(x)
bc.d[prune==0] <- 0

#Step 3. Reduce to percolation threshold
thresh <- percThreshold(bc.d)$threshold
pruned.net <- bc.d
pruned.net[bc.d<thresh] <- 0

if (plot.net){
  v.cex <- apply(x,2,sum) #scaling node size by the log of species frequencies
  v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*scalar)+min.vsize
  gplot(abs(pruned.net),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.lwd=(abs(pruned.net)),vertex.cex=v.cex,vertex.col='grey')
}

return(pruned.net)

}

###Araujo network statistics

#species strength
coStrength <- function(x='dependency network',direction='in'){
  if (direction=='in'){
    return(apply(x,2,sum))    
  }else{
    return(apply(x,1,sum))
  }
}

#symmetry
coSym <- function(x='dependency network',zero.na=TRUE){
  out <- x * 0
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (max(c(x[i,j],x[j,i]))==0){}else{
        out[i,j] <- abs(x[i,j]-x[j,i])/max(c(x[i,j],x[j,i]))
      }
    }
  }
  return(out)
}


###Ecological

###Calculate Stone and Roberts C-score
cscore <- function(x,cu.mat=FALSE){
  x[x!=0] <- 1 #force binary
  cu <- matrix(0,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      ri <- sum(x[,i])
      rj <- sum(x[,j])
      S <- x[,i]*0
      S[x[,i]==1&x[,j]==1] <- 1
      S <- sum(S)
      cu[i,j] <- (ri-S)*(rj-S)
    }
  }
  if (cu.mat){return(cu)}else{return(mean(cu))}
}



nullCom <- function(com,method='r1',nits=5000,burn=500,thin=10){
                                        #force binary
  com[com!=0] <- 1
  for (i in 1:burn){post.burn <- commsimulator(x=com,method=method,thin=thin)}
  out <- list()
  for (i in 1:nits){out[[i]] <- commsimulator(x=post.burn,method=method,thin=thin)}
  return(out)
}


###Evolutionary

##Effective community diversity calculator

##Mass action community genetics simulator
