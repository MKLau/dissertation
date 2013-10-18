###Arthropod Co-occurrence Patterns
###Gina and Art's datasets

library(xlsx)
library(vegan)
library(ecodist)
library(rgl)
library(mvabund)
library(sna)
source('~/cor_nets/CorNets.R')

###Gina's data
                                        #using 2003 data only
gina <- read.xlsx('~/data/2000-2003 Garden Data from Gina.xls',4)
                                        #remove backcrosses
gina <- gina[substr(gina[,1],1,2) != 'bc',]
                                        #make NAs zeros
gina[is.na(gina)] <- 0
                                        #mvabund
y <- mvabund(gina[,-1:-2])
y <- y[,apply(y,2,sum)!=0]
x <- factor(substr(gina[,1],1,2))
gina.mglm <- manyglm(y~x)
plot(gina.mglm)
gina.maov <- anova(gina.mglm)
p.adjust(gina.maov$uni.p[2,],method='fdr')

#networks
gina.anal <- co.nets(cbind(substr(gina[,1],1,2),gina[,-1:-2]),min.dens=3)
net <- gina.anal$cp.net
v.cex <- (degree(net)/max(degree(net))+1)^2
gplot(net,gmode='graph',displaylabels=TRUE,label.cex=0.5,vertex.sides=50,vertex.cex=v.cex,edge.col='grey',edge.lwd=0.65,vertex.col='darkorange',vertex.border='lightblue')
gina.anal$spp.names[v.cex==max(v.cex)]
colnames(net)[v.cex==max(v.cex)]
gina.anal$spp.names[c(26,55,71,30,48,77,9)+3]
#modules
detach(package:sna)
library(igraph)
main.net. <- main.net <- net[apply(net,1,sum)!=0,apply(net,2,sum)!=0]
main.net.[main.net.!=0] <- 1
net.igraph <- graph.adjacency(main.net.)
gina.clust <- clusters(net.igraph,'strong')$membership+1
detach(package:igraph)
library(sna)
v.col <- heat.colors(max(gina.clust))[gina.clust]
v.col <- grey(gina.clust/max(gina.clust))
gplot(main.net,gmode='graph',displaylabels=TRUE,label.cex=0.5,vertex.sides=50,vertex.cex=v.cex,edge.col='grey',edge.lwd=0.65,vertex.col=v.col,vertex.border='lightblue')

###Art's Data
  art <- list()
art[[1]] <- read.csv('~/data/onc_arth_04.csv')
art[[2]] <- read.csv('~/data/onc_arth_05.csv')
art[[3]] <- read.csv('~/data/onc_arth_06.csv')
x <- art[[2]]
                                        #analysis
art.anal <- lapply(art,function(x) co.nets(x,min.dens=5))
names(art.anal) <- 2004:2006
                                        #plots
par(mfrow=c(1,3))
for (i in 1:length(art)){
  net <- art.anal[[i]]$cp.net
  v.cex <- (degree(net)/max(degree(net))+1)^2
  e.col <- net
  e.col[e.col>0.5] <- 'black'
  e.col[e.col==0.5] <- 'grey'
  e.col[e.col<0.5] <- 'red'
  gplot(net,gmode='graph',vertex.sides=50,displaylabels=TRUE,label.cex=0.75,vertex.col='grey',vertex.cex=v.cex,edge.col=e.col)
  title(main=names(art.anal)[i])
}
                                        #correlation between degree and frequency
my.degree <- function(x){
  diag(x) <- 0
  x[x!=0] <- 1
  x <- apply(x,1,sum)
  return(x)
}
par(mfrow=c(1,3))
plot(my.degree(art.anal[[1]]$dist.net)~apply(art.anal[[1]]$x,2,sum))
plot(my.degree(art.anal[[2]]$dist.net)~apply(art.anal[[2]]$x,2,sum))
plot(my.degree(art.anal[[3]]$dist.net)~apply(art.anal[[3]]$x,2,sum))
                                        #
summary(glm(my.degree(art.anal[[1]]$dist.net)~apply(art.anal[[1]]$x,2,sum),family='poisson'))
summary(glm(my.degree(art.anal[[2]]$dist.net)~apply(art.anal[[2]]$x,2,sum),family='poisson'))
summary(glm(my.degree(art.anal[[3]]$dist.net)~apply(art.anal[[3]]$x,2,sum),family='poisson'))
                                        #
par(mfrow=c(2,2))
plot(glm(my.degree(art.anal[[1]]$dist.net)~apply(art.anal[[1]]$x,2,sum),family='poisson'))
plot(glm(my.degree(art.anal[[2]]$dist.net)~apply(art.anal[[2]]$x,2,sum),family='poisson'))
plot(glm(my.degree(art.anal[[3]]$dist.net)~apply(art.anal[[3]]$x,2,sum),family='poisson'))

                                        #
n <- c(47:60)
art.anal[[2]]$spp.names[(n+3)]
data.frame(n,art.anal[[2]]$spp.names[(n+3)])
n <- 7
art.anal[[2]]$spp.names[(n+3)]
(1:length(art.anal[[2]]$spp.names))[art.anal[[2]]$spp.names=='PB.Gall']+3
art.anal[[2]]$spp.names[(46+3)]
                                        #mvabund
art.maov <- list()
min.spp <- 5
                                        #2004
x <- art[[1]][,1]
y <- mvabund(art[[1]][,-1])
y <- y[,apply(y,2,sum) > min.spp]
mva <- manyglm(y~x)
par(mfrow=c(1,3))
plot.manyglm(mva)
meanvar.plot(y~x)
maov <- anova(mva,p.uni='unadjusted')
names(maov)
art.maov[[1]] <- maov
                                        #2005
x <- art[[2]][,1]
y <- mvabund(art[[2]][,-1])
y <- y[,apply(y,2,sum) > min.spp]
mva <- manyglm(y~x)
par(mfrow=c(1,3))
plot.manyglm(mva)
meanvar.plot(y~x)
maov <- anova(mva,p.uni='unadjusted')
names(maov)
maov$uni.p
art.maov[[2]] <- maov
                                        #2006
x <- art[[3]][,1]
y <- mvabund(art[[3]][,-1])
y <- y[,apply(y,2,sum) > min.spp]
mva <- manyglm(y~x)
par(mfrow=c(1,3))
plot.manyglm(mva)
meanvar.plot(y~x)
maov <- anova(mva,p.uni='unadjusted')
art.maov[[3]] <- maov
                                        #compare years
                                        #FDR
p.adjust(art.maov[[1]]$uni.p[2,],method='fdr')[p.adjust(art.maov[[1]]$uni.p[2,],method='fdr')<=0.10]
p.adjust(art.maov[[2]]$uni.p[2,],method='fdr')[p.adjust(art.maov[[2]]$uni.p[2,],method='fdr')<=0.10]
p.adjust(art.maov[[3]]$uni.p[2,],method='fdr')[p.adjust(art.maov[[3]]$uni.p[2,],method='fdr')<=0.10]
                                        #unadjusted
art.maov[[1]]$uni.p[2,][art.maov[[1]]$uni.p[2,]<=0.05]
art.maov[[2]]$uni.p[2,][art.maov[[2]]$uni.p[2,]<=0.05]
art.maov[[3]]$uni.p[2,][art.maov[[3]]$uni.p[2,]<=0.05]
                                        #
n <- 19
art.anal[[2]]$spp.names[art.anal[[2]]$spp.names==art.anal[[2]]$spp.names[n+9]]
n <- 3
(1:length(art.anal[[n]]$spp.names))[art.anal[[n]]$spp.names=='PB.Gall']-9
library(gplots)
mu <- tapply(art[[1]][,colnames(art[[1]])=='PB.Gall'],art[[1]][,1],mean)
se <- tapply(art[[1]][,colnames(art[[1]])=='PB.Gall'],art[[1]][,1],function(x) sd(x)/length(x))
mu <- rbind(mu,tapply(art[[2]][,colnames(art[[2]])=='PB.Gall'],art[[2]][,1],mean))
se <- rbind(se,tapply(art[[2]][,colnames(art[[2]])=='PB.Gall'],art[[2]][,1],function(x) sd(x)/length(x)))
mu <- rbind(mu,tapply(art[[3]][,colnames(art[[3]])=='PB.Gall'],art[[3]][,1],mean))
se <- rbind(se,tapply(art[[3]][,colnames(art[[3]])=='PB.Gall'],art[[3]][,1],function(x) sd(x)/length(x)))
barplot2(mu,beside=TRUE,las=2,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se)

