###LCO - Analysis of the wild stand co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##Site = Uintah, UT
##Study area = 225 * 463 = 104,175 m2 = 0.104175 km2



###lichen data
x <- read.csv('~/projects/dissertation/projects/lichen_coo/data/lco_Apr2012.csv')
                                        #remove notes
x <- x[,colnames(x)!='NOTES.']
                                        #
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
                                        #condense species
lec.spp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x,lec=lec.spp,phy=phy.spp)
x <- x[,-c(6,8,10,18,13,14,15,16)]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
x.q <- split(x,quads)
                                        #data from lamit
env <- read.csv('~/projects/dissertation/projects/lichen_coo/data/Uinta2012_all_data_from_Lamit.csv')
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
x.split <- paste(x$tree,x$quadrat,sep='_')
env.split <- paste(env$Tree.ID,env$Quad.Loc)
x.split <- as.character(x$tree)
env.split <- as.character(env$Tree.ID)
prb <- tapply(env$Pct.Roughness,env.split,mean) #percent rough bark
                                        #environmental data from lamit
wenv <- read.csv('~/projects/dissertation/projects/lichen_coo/data/uintah_trees/uintah_2012_lco_trees.csv')
wenv$tree <- sub('\\.0','\\.',wenv$tree) #remove leading zeros
wenv$tree[wenv$tree=='gnu.468'] <- "gnu.68.ramet"
                                        #match to lichen dataset
wenv <- wenv[wenv$tree%in%names(prb),]
wenv <- wenv[match(names(prb),wenv$tree),]
                                        #lat and long
lat <- wenv$lat..N.deg_decimal.
lat <- t(sapply(lat,function(x) as.numeric(unlist(strsplit(as.character(x),split='_')))))
lat <- lat[,1] + lat[,2]/60
lon <- wenv$lon
lon <- t(sapply(lon,function(x) as.numeric(unlist(strsplit(as.character(x),split='_')))))
lon <- lon[,1] + lon[,2]/60
wenv$lat <- lat
wenv$lon <- lon
                                        #calculate geographic distance
library(fossil)
geo.dist <- earth.dist(data.frame(lon,lat))
range(geo.dist * 1000)

###Tree id definitely influences the communities
###Accounting for tree autocorrelation of observations
com <- split(x,paste(x$tree,x$quadrat,sep='_'))
com <- lapply(com,function(x) apply(x[,-1:-4],2,sum))
com <- do.call(rbind,com)
library(vegan)
com. <- cbind(com,ds=rep(1,nrow(com)))
tid <- sapply(rownames(com.),function(x) strsplit(x,split='_')[[1]][1])
adonis(com.~tid)
###Height doesn't influence communities
ht <- sapply(rownames(com.),function(x) strsplit(x,split='_')[[1]][2])
adonis(com.~ht)

##Species accumulation curves
sac.com <- list()
for (i in 1:length(unique(env$Tree.ID))){
  sac.com[[i]] <- apply(com.[env$Tree.ID==unique(env$Tree.ID)[i],],2,sum)
}
sac.com <- do.call(rbind,sac.com)
plot(specaccum(sac.com),add=FALSE,col=2)

###Co-occurrence C-score

##Whole Stand
source('../src/seenetR.R')
library(sna)
scn <- CoNetwork(x[,((1:ncol(x))[colnames(x)=='xgal']):((1:ncol(x))[colnames(x)=='phy'])])

##Bipartite network (trees and species)
bpn.l <- split(x,as.character(x$tree))
bpn.l <- lapply(bpn.l,function(x) x[,-1:-4])
bpn <- do.call(rbind,lapply(bpn.l,function(x) apply(x,2,sum)))
bpn.sort <- bpn[order(apply(bpn,1,function(x) sum(sign(x))),decreasing=TRUE),order(apply(bpn,2,function(x) sum(sign(x))),decreasing=TRUE)]
plotweb(bpn.sort,method='normal')
nest.r0 <- oecosimu(bpn,nestedtemp,method='r0',nsimul=1000,burnin=100)
nest.r1 <- oecosimu(bpn,nestedtemp,method='r1',nsimul=1000,burnin=100)

##Centrality analysis
detach(package:igraph)
scn.centrality <- degree(scn)
scn.abund <- apply(x[,((1:ncol(x))[colnames(x)=='xgal']):((1:ncol(x))[colnames(x)=='phy'])],2,sum)
                                        #
par(mfrow=c(1,2))
plot(log(scn.centrality+1)~scn.abund)
abline(lm(log(scn.centrality+1)~scn.abund))
summary(lm(log(scn.centrality+1)~scn.abund))
plot(scn.centrality[scn.centrality!=0]~scn.abund[scn.centrality!=0])
abline(lm(scn.centrality[scn.centrality!=0]~scn.abund[scn.centrality!=0]))
summary(lm(scn.centrality[scn.centrality!=0]~scn.abund[scn.centrality!=0]))
                                        #
##Condensed stand level
x.t <- split(x,x.split)
x.t <- split(x,x.split)
x.t <- lapply(x.t,function(x) x[,-1:-4])
com <- do.call(rbind,lapply(x.t,function(x) apply(x,2,sum)))
cond <- com
cond[cond!=0] <- 1
cond.ncs <- nullCom(com)
cond.ncs <- unlist(lapply(cond.ncs,cscore))
cond.ocs <- cscore(com)
cond.ses <- (cond.ocs-mean(cond.ncs))/sd(cond.ncs)
cond.p <- length(cond.ncs[cond.ncs<=cond.ocs])/length(cond.ncs)
cond.net <- dep.net(cond)
adonis(com~prb,permutations=10000)
com.rel <- apply(com,2,function(x) x/max(x))
adonis(com.rel~prb,permutations=10000)
                                        #network stats from Araujo
net.sim <- coSym(dep.net(x[,((1:ncol(x))[colnames(x)=='xgal']):((1:ncol(x))[colnames(x)=='phy'])]))
mean(net.sim[upper.tri(net.sim)])

                                        #max network plot
max.order <- order(prb,decreasing=TRUE)
net <- scn
max.com <- com[max.order,]
rownames(max.com) <- 1:nrow(max.com)
max.loc <- function(x){
  x[x!=max(x)] <- 0
  x[x==max(x)] <- max(x)
  return(x)
}
max.t <- apply(max.com,2,max.loc)
max.pos <- rep(3,ncol(max.t))
max.xy <- array(NA,dim=c(ncol(max.t),2))
rownames(max.xy) <- colnames(max.t)
                                        #plotting
par(mfrow=c(1,1))
plot(c(0,max(max.t))~c(0,nrow(max.t)+0.5),type='n',ylab='Maximum Abundance',xaxt='n',xlab='Tree (Ranked by Roughness)',font.lab=2)
axis(1,at=seq(0.5,nrow(max.t),length=nrow(max.t)),label=rownames(max.t),las=1)
for (i in 1:ncol(max.t)){
  x <- seq(0.5,nrow(max.t),length=nrow(max.t))[max.t[,i]!=0]
  if (i==6){x <- x + 0.2}
  y <- max.t[max.t[,i]!=0,i]
  max.xy[i,] <- c(x,y)
  points(x,y,pch=19)
  text(x,y,colnames(max.t)[i],pos=max.pos[i])
}
                                        #networking
for (i in 1:nrow(net)){
  for (j in 1:ncol(net)){
    if (net[i,j]!=0){
      lines(c(max.xy[rownames(max.xy)==rownames(net)[i],1],max.xy[rownames(max.xy)==colnames(net)[j],1]),
            c(max.xy[rownames(max.xy)==rownames(net)[i],2],max.xy[rownames(max.xy)==colnames(net)[j],2]),
            col='grey')
    }
  }
}

##Tree scale
cs <- unlist(lapply(x.t,cscore))
tid <- unlist(sapply(names(cs),function(x) strsplit(x,split='_')[[1]][1]))
                                        #co-occurrence nets
source('../src/seenetR.R')
dn.t <- lapply(x.t,CoNetwork) #Co-occurrence networks
net.d <- matrix(0,nrow=length(dn.t),ncol=length(dn.t))
rownames(net.d) <- colnames(net.d) <- names(dn.t)
for (i in 1:nrow(net.d)){
  for (j in 1:ncol(net.d)){
    net.d[i,j] <- sum(abs(dn.t[[i]]-dn.t[[j]])^2)
  }
}
net.d <- as.dist(net.d)
                                        #averaged network structure
avg.dnet <- dn.t[[1]]
for (i in 2:length(dn.t)){
  avg.dnet <- avg.dnet + dn.t[[i]]
}
avg.dnet <- avg.dnet / length(avg.dnet)

par(mfrow=c(1,2))
v.cex <- apply(x[-1:-4],2,sum) #scaling node size by the log of species frequencies
v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*3)+0.1
coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
               edge.lwd=(abs(scn)),vertex.cex=v.cex,vertex.col='grey')
gplot(abs(avg.dnet),displaylabels=TRUE,gmode='graph',pad=1.5,
      edge.lwd=(abs(avg.dnet)),vertex.cex=v.cex,vertex.col='grey',
      coord=coord)

                                        #wild ses scores
ses.t <- dget('../data/wild_tree_ses.Rdata')
ses.t[is.na(ses.t)] <- 0
ses.p <- dget('../data/wild_tree_pval.Rdata')
ses.ncs <- dget('../data/wild_tree_ncs.Rdata')
ses.ocs <- dget('../data/wild_tree_ocs.Rdata')
names(ses.t) <- names(ses.ocs)
                                        #centrality
detach(package:igraph)
cen.t <- unlist(lapply(dn.t,function(x) centralization(x,FUN='degree')))
                                        #degree
deg.t <- unlist(lapply(dn.t,function(x) sum(sign(abs(x)))))
                                        #distance
mantel(geo.dist,net.d,permutations=5000) #environmental distance
adonis(net.d~prb)
adonis(net.d~deg.t)
adonis(net.d~cen.t)
summary(lm(cen.t~deg.t))
summary(lm(deg.t~prb))
summary(lm(cen.t~prb))
summary(lm(ses.t~deg.t))
summary(lm(ses.t~prb))
                                        #SES values
###SES roughness, age, abundance, richness
height <- as.character(sapply(rownames(com),function(x) strsplit(x,split='_')[[1]][2]))
com.4555 <- com[height=='n45.55',]
com.8090 <- com[height=='n80.90',]
com.t <- com.4555+com.8090
rownames(com.t) <- as.character(unlist(sapply(rownames(com.4555),function(x) strsplit(x,split='_')[[1]][1])))
total.abundance <- apply(com.t,1,sum)
species.richness <- apply(sign(com.t),1,sum)
                                        #composition by rough
adonis((com.t)~prb,permutation=5000)
adonis((com.t)~age,permutation=5000)
mantel(vegdist(com.t)~geo.dist)
mantel(vegdist(com.t)~net.d+geo.dist)
                                        #
plot(ses.t~prb,xlab='Mean Percent Roughness',ylab='SES',pch=19,
     cex=1.5,font.lab=2,cex.lab=1.25,font.axis=2)
my.line(prb,ses.t,lwd=3)
summary(lm(ses.t~prb))
                                        #
plot((ses.t~total.abundance),xlab='Total Abundance',ylab='SES',pch=19,font.lab=2)
abline(lm(ses.t~total.abundance))
summary(lm(ses.t~total.abundance))
                                        #
plot((ses.t~species.richness),xlab='Species Richness',ylab='SES',pch=19,font.lab=2)
abline(lm(ses.t~species.richness))
summary(lm(ses.t~species.richness))
                                        #Ageage <- read.csv('~/projects/dissertation/projects/lichen_coo/data/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv')
                                        #
age <- read.csv('~/projects/dissertation/projects/lichen_coo/data/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv')
dbh <- age$DBH.cm_01
age.final <- age$AgeFinal.U
age <- data.frame(tree.id=age[,1],age.final=age$AgeFinal.U)
age[,1] <- tolower(age[,1])
age[,1] <- sub('_','\\.',age[,1])
age[,1] <- sub('-','\\.',age[,1])
age[,1] <- sub('\\?','',age[,1])
age[,1] <- sub('\\.0','\\.',age[,1])
age[age[,1]=='gnu.85.1ftaway',1] <- 'gnu.85'
                                        #
                                        #predict age
gnu19.dbh <- dbh[age$tree.id=='gnu.19']
new <- data.frame(dbh=seq(min(dbh),max(dbh),by=0.1))
age.final <- na.omit(age.final)
pred.age <- predict(lm(age.final~dbh,data=age),new)
plot(pred.age~new[,1])
gnu19.age <- as.numeric(pred.age[new[,1]==gnu19.dbh])
                                        #
tree.age <- numeric(length(ses.t))
tree.age <- age[match(names(ses.t),age[,1]),2]
tree.age[is.na(tree.age)] <- gnu19.age
names(tree.age) <- age[match(names(ses.t),age[,1]),1]
                                        #
par(mfrow=c(1,2))
plot(prb~tree.age,xlab='Tree Age (years)',ylab='Percent Bark Roughness',font.lab=2,pch=19)
abline(lm(prb~tree.age))
summary(lm(prb~tree.age))
plot(ses.t~tree.age)
abline(lm(ses.t~tree.age))
summary(lm(ses.t~tree.age))
                                        #
library(sem)
###tree.age -> prb -> ses
model.sem <- specifyModel(file='~/projects/dissertation/projects/lichen_coo/src/sem_model.txt')
   sem.data <- na.omit(cbind(tree.age,prb,ses=ses.t))
                                        #transforms
###
                                        #sem fitting and analysis
   Sigma <- var(sem.data)
   sem.fit <- sem(model.sem,S=Sigma,N=nrow(sem.data))
   summary(sem.fit)

                                        #effects(sem.fit)
                                        #hist(residuals(sem.fit))
                                        #stdCoef(sem.fit)
pathDiagram(sem.fit,file='~/projects/dissertation/projects/lichen_coo/results/sem_coo',
               edge.labels='values',standardize=TRUE
               ,ignore.double=TRUE,size=c(12,12),edge.font=c("Arial", 10),
               graphics.fmt='png') #export to graphviz

###Tree characters predict network structure (deg.t) and co-occurrence patterns (ses.zp)
                                        #roughness (prb)
summary(lm(deg.t~prb))
summary(lm(ses.t~prb))
                                        #cover (canopy.n)/96
cover <- (wenv$canopy.n/96)*100
summary(lm(cover~tree.age))
summary(lm(deg.t~cover))
summary(lm(ses.t~cover))
                                        #height
summary(lm(ses.t~wenv$hieght.m))

###################################
###Final analyses
###################################
                                        #stand network
par(mfow=c(1,1))
gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
      edge.col=e.col,edge.lwd=abs(scn),
      vertex.cex=v.cex,vertex.col='lightblue')
                                        #tree scale
                                        #affect of PRB on SES
summary(lm(ses.t~prb))
shapiro.test(residuals(summary(lm(ses.t~prb))))
                                        #affect of age on ses
summary(lm(ses.t~tree.age))
shapiro.test(residuals(summary(lm(ses.t~tree.age))))
                                        #mantels
                                        #microsite
                                        # geographic distance
                                        # elevation
mantel(dist(ses.t)~dist(prb)+dist(tree.age)+geo.dist)
mantel(dist(ses.t)~dist(tree.age)+geo.dist+dist(prb))
mantel(dist(prb)~dist(tree.age)+geo.dist)
mantel(dist(prb)~geo.dist)
mantel(dist(tree.age)~geo.dist)
                                        #SEM for reduced model
library(sem)
###tree.age -> prb -> ses
model.sem <- specifyModel(file='~/projects/dissertation/projects/lichen_coo/src/sem_ses.txt')
   sem.data <- cbind(tree.age,prb,ses,cover)
                                        #transforms
###
                                        #sem fitting and analysis
   Sigma <- var(sem.data)
   sem.fit <- sem(model.sem,S=Sigma,N=nrow(sem.data))
   summary(sem.fit)
   modIndices(sem.fit)
                                        #effects(sem.fit)
                                        #hist(residuals(sem.fit))
                                        #stdCoef(sem.fit)
pathDiagram(sem.fit,file='~/projects/dissertation/projects/lichen_coo/results/sem_coo',
               edge.labels='values',standardize=TRUE
               ,ignore.double=TRUE,size=c(12,12),edge.font=c("Arial", 10),
               graphics.fmt='png') #export to graphviz
