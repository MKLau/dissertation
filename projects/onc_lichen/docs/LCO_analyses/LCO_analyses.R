###################################################
### chunk number 1: 
###################################################
#line 37 "LCO_analyses.Rnw"
require(RLRsim)
require(lme4)
require(sna)
#require(network)
#require(igraph)                                      
require(vegan)
require(xtable)
require(gplots)
source('source/CorNets.R')
source('source/pairs.R')
se = function(x){sd(x)/sqrt(length(x))}

binary = function(x,bin=0){x[x>bin]=1;return(x)}

richness = function(x){
  x = apply(x,2,sum)
  x[x != 0] = 1
  sum(x)
}

bplot = function(x,y,ylab=''){
 mu. = tapply(y,x,mean)
 se. = tapply(y,x,se)
 barplot2(mu.,plot.ci=TRUE,ci.u=mu.+se.,ci.l=mu.-se.,las=2,ylab=ylab)
}

which.garden = function(x){
 x = unlist(strsplit(as.character(x),split='.',fixed=TRUE))[1]
 x = unlist(strsplit(x,split=''))
 if (length(x) < 3){x = 'ONC'}else {x = 'PIT'}
 return(x)
}

lco2quad = function(x){
 com = x[,-1:-6]
 com = apply(com,1,sum)
 x = x[,5:6]
 out = array(0,dim=c(max(x),max(x)))
 for (i in (1:nrow(x))){
  out[x[i,1],x[i,2]] = com[i]
 }
 return(out)
}

exp.m <- function(mat, n){
  ## Alberto Monteiro
  ## https://stat.ethz.ch/pipermail/r-help/2007-May/131330.html
  ## test if mat is a square matrix
  ## treat n < 0 and n = 0 -- this is left as an exercise
  ## trap non-integer n and return an error
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat %*% mat
    n <- n / 2
  }
  return(result)
}

proliferate = function(g = 'graph',lim = 10000){
  n = numeric()
  for (i in 1:lim){
    n[i] = sum(exp.m(g,i))
    if (n[i] == Inf){break}
  }
n = n[-i]
return(n)
}



###################################################
### chunk number 2: 
###################################################
#line 116 "LCO_analyses.Rnw"
  lco = read.csv('data/LCO_data_ONC_PIT.csv')


###################################################
### chunk number 3: 
###################################################
#line 121 "LCO_analyses.Rnw"
  xtable(summary(lco[,1:6]),table.placement='tbp')


###################################################
### chunk number 4: 
###################################################
#line 125 "LCO_analyses.Rnw"
  xtable(summary(lco[,7:10]),table.placement='tbp')


###################################################
### chunk number 5: 
###################################################
#line 129 "LCO_analyses.Rnw"
  xtable(summary(lco[,13:ncol(lco)]),table.placement='tbp')


###################################################
### chunk number 6: 
###################################################
#line 134 "LCO_analyses.Rnw"
  lco = lco[lco$Quadrat == 'n45.55',]  #remove upper quadrats
  lco.l = list()
  tree.names = as.character(unique(lco$Tree))


###################################################
### chunk number 7: 
###################################################
#line 140 "LCO_analyses.Rnw"
  for (i in (1:length(tree.names))){
    lco.l[[i]] = lco[lco$Tree == tree.names[i],]
  }
  names(lco.l) = tree.names
  any(table(tree.names) != 1)  #check for any duplicate trees


###################################################
### chunk number 8: 
###################################################
#line 149 "LCO_analyses.Rnw"

  garden = sapply(tree.names,which.garden)
  geno = unlist(lapply(lco.l,function(x) as.character(x$Geno[1])))
  com.l = lapply(lco.l,function(x) x[,(1:ncol(x))[colnames(x) == 'Xgal']:ncol(x)])



###################################################
### chunk number 9: 
###################################################
#line 157 "LCO_analyses.Rnw"
geno.table = tapply(factor(geno),garden,table)
geno.table = rbind(geno.table$'ONC',geno.table$'PIT')
rownames(geno.table) = c('ONC','PIT')
geno.table = geno.table[,order(apply(geno.table,2,sum),decreasing=TRUE)]
xtable(geno.table)
barplot(geno.table,beside=FALSE,las=2,ylab='Number of Trees',col=c(1,0))
legend('topright',c('ONC','PIT'),fill=c(1,0),bty='n')


###################################################
### chunk number 10: 
###################################################
#line 169 "LCO_analyses.Rnw"

 quad.l = lapply(lco.l,lco2quad)



###################################################
### chunk number 11: 
###################################################
#line 179 "LCO_analyses.Rnw"
zero.na = function(x){x[is.na(x)] = 0;return(x)}
cor.l = lapply(com.l,kendall.pairs, alpha = 0.05,p.adj=TRUE)
cor.l = lapply(cor.l,zero.na)


###################################################
### chunk number 12: 
###################################################
#line 187 "LCO_analyses.Rnw"
  A = lapply(com.l,function(x) apply(x,2,sum))
  A.= lapply(A,binary)
  R = unlist(lapply(com.l,richness))
  n = unlist(lapply(cor.l,function(x) length(apply(x,1,sum)[abs(apply(x,1,sum)) > 0]))) #number of nodes
  L = unlist(lapply(cor.l,function(x) length(x[x != 0])/2)) #number of edges
  C = L/(n*(n-1)); C[is.na(C)] = 0 #connectance
  abs.cor.l = lapply(cor.l,abs)
  C.f = unlist(lapply(abs.cor.l,function(x) centralization(x,degree,mode='graph')))  #Freeman's centralization




###################################################
### chunk number 13: 
###################################################
#line 242 "LCO_analyses.Rnw"
input.pairs = cbind(R,n,L,C,C.f)
pairs.onc = input.pairs[garden == 'ONC',]
pairs.pit = input.pairs[garden == 'PIT',]



###################################################
### chunk number 14: 
###################################################
#line 249 "LCO_analyses.Rnw"
pairs(input.pairs,upper.panel=panel.lm,lower.panel=panel.cor)



###################################################
### chunk number 15: 
###################################################
#line 254 "LCO_analyses.Rnw"
pairs(pairs.onc,upper.panel=panel.lm,lower.panel=panel.cor)



###################################################
### chunk number 16: 
###################################################
#line 259 "LCO_analyses.Rnw"
pairs(pairs.pit,upper.panel=panel.lm,lower.panel=panel.cor)



###################################################
### chunk number 17: 
###################################################
#line 268 "LCO_analyses.Rnw"
##ONC Analyses
 geno.onc = geno[garden == 'ONC']
 n.onc = n[garden == 'ONC']
 log.n.onc = log(n.onc+0.5)
 L.onc = L[garden == 'ONC']
 log.L.onc = log(L.onc+0.5)
 C.onc = C[garden == 'ONC']
 Cf.onc = C.f[garden == 'ONC']
 asin.sqrt.CF.onc = asin(sqrt(Cf.onc))
                                        #Size
 exactRLRT(lmer(n.onc~1|geno.onc)) 
 exactRLRT(lmer(log.n.onc~1|geno.onc)) 
                                        #Number of Edges
 exactRLRT(lmer(L.onc~1|geno.onc)) 
 exactRLRT(lmer(log.L.onc~1|geno.onc)) 
                                        #Connectance
 exactRLRT(lmer(C.onc~1|geno.onc)) 
                                        #Freeman's Degree Centrality
 exactRLRT(lmer(Cf.onc~1|geno.onc))
 exactRLRT(lmer(asin.sqrt.CF.onc~1|geno.onc))                                      

 geno.pit = geno[garden == 'PIT'] 
 n.pit = n[garden == 'PIT']
 log.n.pit = log(n.pit+0.5)
 L.pit = L[garden == 'PIT']
 log.L.pit = log(L.pit+0.5)
 C.pit = C[garden == 'PIT']
 Cf.pit = C.f[garden == 'PIT']
 asin.sqrt.CF.pit = asin(sqrt(Cf.pit))
                                        #Size
 exactRLRT(lmer(n.pit~1|geno.pit))
 exactRLRT(lmer(log.n.pit~1|geno.pit))
                                        #Number of Edges
 exactRLRT(lmer(L.pit~1|geno.pit))
 exactRLRT(lmer(log.L.pit~1|geno.pit))
                                        #Connectance
 exactRLRT(lmer(C.pit~1|geno.pit))
                                        #Freeman's Degree of Centralization
 exactRLRT(lmer(Cf.pit~1|geno.pit))
 exactRLRT(lmer(asin.sqrt.CF.pit~1|geno.pit))


###################################################
### chunk number 18: 
###################################################
#line 333 "LCO_analyses.Rnw"
  gplot(G.graph,gmode='graph',displaylabels=TRUE,mode='circle',edge.lwd=G.graph*200)


###################################################
### chunk number 19: 
###################################################
#line 340 "LCO_analyses.Rnw"
cor.l. = cor.l[garden == 'ONC']
geno.A = A[garden == 'ONC']
geno.A.= A.[garden == 'ONC']
geno.net = list()
for (i in (1:length(unique(geno.onc)))){
  x = cor.l.[geno.onc == unique(geno.onc)[i]]
  geno.net[[i]] = x[[1]]
  for (j in (2:length(x))){
    geno.net[[i]] = geno.net[[i]] + x[[j]]
  }
  geno.net[[i]] = geno.net[[i]] / length(x)
}
names(geno.net) = unique(geno)


###################################################
### chunk number 20: 
###################################################
#line 356 "LCO_analyses.Rnw"
par(mfrow=c(4,4))
for (i in (1:length(geno.net))){
gplot(abs(geno.net[[i]]),gmode='graph',displaylabels=FALSE,mode='circle',main=names(geno.net)[i],cex.main=2,edge.lwd=abs(geno.net[[i]])*100,vertex.col=geno.A.[[i]],edge.col='grey',vertex.sides=100,vertex.cex=(geno.A[[i]]/max(geno.A[[i]])+1)^2,vertex.lty=geno.A.[[i]])
}


###################################################
### chunk number 21: 
###################################################
#line 365 "LCO_analyses.Rnw"
 par(mfrow=c(2,2))
 bplot(geno.onc,n.onc,'Size')
 bplot(geno.onc,L.onc,'Number of Edges')
 bplot(geno.onc,C.onc,'Connectance')
 bplot(geno.onc,Cf.onc,'Centralization')


###################################################
### chunk number 22: 
###################################################
#line 376 "LCO_analyses.Rnw"
 par(mfrow=c(2,2))
 bplot(geno.pit,n.pit,'Size')
 bplot(geno.pit,L.pit,'Number of Edges')
 bplot(geno.pit,C.pit,'Connectance')
 bplot(geno.pit,Cf.pit,'Centralization')


