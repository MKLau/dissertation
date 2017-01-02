###Sunset Lichen Networks

library(vegan)
library(bipartite)
                                        #
key <- read.csv('../data/2012jul30/Key.csv')
x <- read.csv('/Users/Aeolus/projects/dissertation/projects/sunset_lichen/data/2012jul30/spp_env_combined.csv')
                                        #remove dead
x <- x[x$Live.Dead == 1,]
com <- x[,((1:ncol(x))[colnames(x) == 'Acacon']):((1:ncol(x))[colnames(x) == 'Xanele'])]
env <- x[,1:12]
                                        #remove N and S light
env <- env[,-10:-11]
                                        #fix colnames
colnames(env) <- sub('\\.\\.','',colnames(env))
colnames(env)

###Composition analysis
com. <- cbind(com,ds=rep(0.01,nrow(com)))
adonis(com.~env$Moth)

###bipartite network
env. <- env[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE),]
bpn <- com[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE),order(apply(com,2,function(x) sum(sign(x))),decreasing=TRUE)]
rownames(bpn) <- paste(env.$Tree.pairs,env.$Moth,sep='_')
                                        #nest.test <- oecosimu(bpn, nestfun="nestedtemp", "r1",nsimul = 999,burnin=50)
nest.test <- dget(file='../results/nest_test.Rdata')
plotweb(bpn,method='normal',text.rot=90,col.low=env.$Moth+1)

###Modularity
bpn.mods <- computeModules(bpn[,apply(bpn,2,function(x) sum(sign(x)))>1])
plotModuleWeb(bpn.mods)

###Araujo Networks
net <- CoNetwork(com[,apply(com,2,sum)!=0])
n.com <- nullCom(com[,apply(com,2,sum)!=0])
n.c <- unlist(lapply(n.com,C.score))
o.c <- C.score(com[,apply(com,2,sum)!=0])
p.val <- length(n.c[n.c<=o.c])/length(n.c)
ses <- (o.c-mean(n.c))/sd(n.c)

v.cex <- apply(com[,apply(com,2,sum)!=0],2,sum) #scaling node size by the log of species frequencies
v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*3)+0.1
e.col <- net
e.col[net>0.5] <- 'red'
e.col[net<0.5] <- 'black'
e.col[net==0.5] <- 'grey'
gplot(abs(net),displaylabels=TRUE,gmode='graph',pad=1.5,
      edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',edge.col=e.col)

###Plot
par(mfcol=c(2,1))
gplot(abs(net),displaylabels=TRUE,gmode='graph',pad=1.5,
      edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',edge.col=e.col)
plotweb(bpn,method='normal',text.rot=90,col.low=env.$Moth+1)

