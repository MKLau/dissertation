###LCO - Analysis of the wild stand co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????

x <- read.csv('~/projects/dissertation/projects/onc_lichen/data/lco_Apr2012.csv')
x <- na.omit(x)
###NOTE: gnu.44 is a fremont
                                        #remove gnu.44
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree
x <- x[x$tree!='ll.6',]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
x.q <- split(x,quads)

###look at abundance and richness patterns
A <- unlist(lapply(x.q,function(x) sum(x[,5:18])))
R <- unlist(lapply(x.q,function(x) sum(sign(apply(x[,5:18],2,sum)))))
plot(R~A)
abline(lm(R~A))
summary(lm(R~A))
height <- unlist(lapply(strsplit(names(x.q),split='\ n'),function(x) x[2]))
summary(aov(A~height))
summary(aov(R~height))

###Composition
library(vegan)
com <- do.call(rbind,lapply(x.q,function(x) apply(x[,5:18],2,sum)))
com. <- cbind(com,ds=rep(1,nrow(com)))
                                        #patterns of height
adonis(com.~height)
com.rel <- apply(com.,2,function(x) if(sum(x)!=0){x/sum(x)}else{x})
adonis(com.rel~height)
                                        #patterns with abundance and richness
adonis(com.~A)
adonis(com.~R)
adonis(com.rel~A)
adonis(com.rel~R)

###Tree identity
tree.id <- sapply(names(x.q),function(x) strsplit(x,split='\ ')[[1]][1])
summary(aov(A~tree.id))
summary(aov(R~tree.id))
adonis(com.~tree.id)
adonis(com.rel~tree.id)

plot(A~as.numeric(factor(tree.id)),pch=19)
plot(R~as.numeric(factor(tree.id)),pch=19)

summary(lm(A[height=='45.55']~R[height=='45.55']))
summary(lm(A[height=='80.90']~R[height=='80.90']))
plot(R~A,col=as.numeric(factor(height)),pch=19)
plot(R~A,col=rainbow(nlevels(factor(tree.id)))[as.numeric(factor(tree.id))],pch=19)

summary(lm(A~R+R/tree.id))

###examine spatial patterns for each species
                                        #reproduce spatial grid for each species
sp.grid <- function(x,sp.col=5){
  x. <- cbind(x[,3:4],x[,sp.col])
  m <- matrix(0,nrow=10,ncol=10)
  for (i in 1:nrow(x.)){
    m[x.[i,1],x.[i,2]] <- x.[i,3]
  }
  return(m)
}
                                        #loops to get the grids
m.q <- list()
m.j <- list()
for (i in 1:length(x.q)){
  for (j in 1:14){
    m.j[[j]] <- sp.grid(x.q[[i]],j+4)
  }
  names(m.j) <- colnames(x.q[[i]])[5:18]
  m.q[[i]] <- m.j
}
                                        #integrated grids
S.q <- lapply(m.q,function(x)x[[1]])
for (i in 1:length(m.q)){
  for (j in 2:length(m.q[[i]])){
    S.q[[i]] <- S.q[[i]] + m.q[[i]][[j]]
  }
}
                                        #plot
image(t(S.q[[i]]),col=heat.colors(length(unique(unlist(S.q)))))

                                        #Null Model Based
cu <- function(x,y){
  S <- sum(x[(x+y)==2])
  return((sum(x)-S)*(sum(y)-S))
}
                                        #
c.score <- function(x,cu.mat=FALSE){
  cu.m <- matrix(0,nrow=ncol(x),ncol=ncol(x))
  rownames(cu.m) <- colnames(cu.m) <- colnames(x)
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      cu.m[i,j] <- cu(x[,i],x[,j])
    }
  }
  if (cu.mat){
    return(cu.m)
  }else{
    return(mean(cu.m))
  }
}
                                        #Cscore analysis
com.q <- lapply(x.q,function(x) x[,5:18])
cs.q <- lapply(com.q,c.score)
cs.q <- unlist(cs.q)
par(mfrow=c(1,2))
plot(cs.q~A)
abline(lm(cs.q~A))
plot(cs.q~R)
abline(lm(cs.q~R))
summary(lm(cs.q~A))
summary(lm(cs.q~R))

summary(lm(cs.q~A+A:tree.id))
summary(lm(cs.q~R+R:tree.id))

par(mfrow=c(1,1))

###Analyses with roughness, age, total cover and height
env <- read.csv('~/projects/dissertation/projects/onc_lichen/data/Uinta2012_all_data_from_Lamit.csv')
env <- env[is.na(env$Pct.Roughness)==FALSE,]
env[,1] <- sub('\\?','',sub('\\.0','\\.',sub('\\_','\\.',sub('\\-','\\.',tolower(as.character(env[,1]))))))
env[env[,1]=='ll.6_(znu.29)',1] <- 'll.6'
env[env[,1]=='gnu.85.1ftaway',1] <- 'gnu.85'
env$Quad.Loc <- as.character(sapply(as.character(env$Quad.Loc),function(x) unlist(strsplit(x,split='_'))[2]))
env$Quad.Loc <- sub('\\-','\\.',env$Quad.Loc)
env$Quad.Loc <- paste('n',env$Quad.Loc,sep='')
                                        #remove southern aspect
env <- env[env$Aspect!='South',]
env.tid <- paste(env$Tree.ID,env$Quad.Loc)
                                        #check that the datasets are compatible
all(names(x.q)%in%env.tid) 
                                        #match observations
all(env.tid[match(names(x.q),env.tid)] == names(x.q))
                                        #delimit to co-occurrence dataset and match
env <- env[match(names(x.q),env.tid),]
                                        #Cscore analyses
all(names(cs.q)==names(x.q))
plot(cs.q~env$Pct.Roughness,font.lab=2,cex.lab=1.25,xlab='Percent Roughness',ylab='C-Score')
abline(lm(cs.q~env$Pct.Roughness))
summary(lm(cs.q~env$Pct.Roughness))
summary(lm(cs.q~env$Pct.Roughness+env$Pct.Roughness:tree.id))

                                        #
plot(cs.q~env$Pct.Total.Cover)
abline(lm(cs.q~env$Pct.Total.Cover))
summary(lm(cs.q~env$Pct.Total.Cover))
summary(lm(cs.q~env$Pct.Total.Cover+env$Pct.Total.Cover:tree.id))
summary(lm(cs.q[height=='45.55']~env$Pct.Total.Cover[height=='45.55']))
summary(lm(cs.q[height=='80.90']~env$Pct.Total.Cover[height=='80.90']))
                                        #isocline for each species

                                        #pairs plots
source('~/projects/dissertation/projects/onc_lichen/docs/LCO_analyses/source/pairs.R')
pairs(cbind(Roughness=env$Pct.Roughness,Total.Cover=env$Pct.Total.Cover,Abundance=A,Richness=R,Cscore=cs.q),upper.panel=panel.lm,lower.panel=panel.cor)
                                        #Look at Araujo co-occurrence networks
library(sna)
source('~/projects/dissertation/projects/onc_lichen/src/araujo_nets.R')

dep.all <- dep.net(all[,5:18])
net.graph <- dep.all[apply(dep.all,1,sum)!=0,apply(dep.all,2,sum)!=0]
gplot(net.graph,displaylabels=TRUE,edge.col='darkgrey',
      gmode='graph',label.cex=1.5,vertex.col='lightblue',
      edge.lwd=(net.graph/max(net.graph)+1)^2.5,uselen=FALSE,
      edge.len=0.025,usearrows=FALSE)
                                        #
net.q <- lapply(x.q,function(x) co.net(x[,5:18]))
unlist(lapply(net.q,function(x) length(x[x!=0])))
                                        #
par(mfrow=c(4,6))
for (i in 1:length(net.q)){
  if (sum(net.q[[i]])==0){
    plot(0)
  }else{
    gplot(net.q[[i]],displaylabels=FALSE)
  }
}

                                        #"pseudo-motifs" = cell species richness = pm
                                        # avg. pseudo.motif size = aps

                                        #SES analysis
                                        #library(pbapply)
                                        #ses.q <- pblapply(com.q,function(x) oecosimu(x[,apply(x,2,sum)!=0], c.score, "swap", burnin=100, thin=10, statistic="evals"))
