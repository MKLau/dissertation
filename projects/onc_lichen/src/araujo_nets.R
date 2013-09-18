###Method for producing network models from co-occurrence data.
###Developed in Araujo et al. 2011 Ecography Using species co-occurrence networks to assess the impacts of climate change
###Coded 18 Sep 2013 MKLau

## E is an environmental lattice with rows = N and cols = M and has size A = N*M.
## A = total number of sites
## Na = number of sites with species a
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

## Three states, i = ab, ii a+b = (a&!b + !a&b), iii = !a!b

## E(i) = APa&b
## E(ii) = APa+b
## E(iii) = AP!a&!b

## Var(i) = APa&b(1-Pa&b)
## Var(ii) = APa+b(1-APa+b)
## Var(iii) = AP!a&!b(1-AP!a&!b)

## Critical interval

## CI.u = E(i) + 2*Var(i)^2 = attractive overlap
## CI.l = E(i) - 2*Var(i)^2 = repulsive overlap

null.prune <- function(a,b,std=TRUE){
                                        #determine number of sites
A <- length(a)
                                        #NULL probabilities
pa <- sum(a)/A
pb <- sum(b)/A
paANDb <- pa * pb
paORb <- pa + pb - 2*pa*pb
pNOTaORb <- 1- pa - pb + pa*pb
                                        #expected frequencies
Ei <- A*paANDb
Eii <- A*paORb
Eiii <- A*pNOTaORb
                                        #variance of frequencies
Vi <- A*paANDb*(1-paANDb)
Vii <- A*paORb*(1-paORb)
Viii <- A*pNOTaORb*(1-pNOTaORb)
                                        #critical threshold (+ = attraction and - = repulstion)
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
dep.net <- function(x='species in cols',zero.na=TRUE,prune=TRUE,diag.zero=TRUE){
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      out[i,j] <- calcDepend(x[,i],x[j])
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
