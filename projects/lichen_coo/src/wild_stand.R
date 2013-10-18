###LCO - Analysis of the wild stand co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##Site = Uintah, UT
##Study area = 225 * 463 = 104,175 m2 = 0.104175 km2

###Tree data
x <- read.csv('~/projects/dissertation/projects/lichen_coo/data/lco_Apr2012.csv')
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
                                        #condense species
lec.spp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x[,-ncol(x)],lec=lec.spp,phy=phy.spp,NOTES=x[,ncol(x)])
x <- x[,-c(6,8,10,18,13,14,15,16)]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
x.q <- split(x,quads)
                                        #environmental data from lamit
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

###Tree id definitely influences the communities
###Accounting for tree autocorrelation of observations
com <- split(x,paste(x$tree,x$quadrat,sep='_'))
com <- lapply(com,function(x) apply(x[,-c(1:4,ncol(x))],2,sum))
com <- do.call(rbind,com)
library(vegan)
com. <- cbind(com,ds=rep(1,nrow(com)))
tid <- sapply(rownames(com.),function(x) strsplit(x,split='_')[[1]][1])
adonis(com.~tid)

###Height doesn't influence communities
ht <- sapply(rownames(com.),function(x) strsplit(x,split='_')[[1]][2])
adonis(com.~ht)

###Co-occurrence C-score

##Whole Stand
source('./seenetR.R')
library(sna)
scn <- dget('~/projects/dissertation/projects/lichen_coo/results/scn.Rdata')
scn <- scn[rownames(scn)!='folio_grey_black',colnames(scn)!='folio_grey_black']
e.col <- sign(scn)
e.col[e.col==1] <- 'grey'
e.col[e.col==-1] <- 'red'
v.cex <- apply(com,2,sum)[colnames(com)!='folio_grey_black']
v.cex <- log(v.cex,10)
v.cex <- v.cex * (1/min(v.cex))
v.cex <- v.cex/2
gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
      edge.col=e.col,edge.lwd=abs(scn),
      vertex.cex=v.cex,vertex.col='lightblue')

##Centrality analysis
scn.centrality <- degree(scn)
scn.abund <- apply(com,2,sum)
scn.rel <- scn.abund/sum(scn.abund)*100
names(scn.rel)[4] <- 'fgb'
                                        #
par(mfrow=c(1,2))
plot(scn.centrality~scn.rel[c(-4)])
abline(lm(scn.centrality~scn.rel[c(-4)]))
summary(lm(scn.centrality~scn.rel[c(-4)]))
                                        #
plot(scn.centrality[-1]~scn.rel[c(-1,-4)])
abline(lm(scn.centrality[-1]~scn.rel[c(-1,-4)]))
summary(lm(scn.centrality[-1]~scn.rel[c(-1,-4)]))
                                        #
par(mfrow=c(1,2))
barplot(scn.centrality[order(scn.centrality,decreasing=TRUE)],
        names=rownames(scn)[order(scn.centrality,decreasing=TRUE)],
        las=2,ylab='Node Level Centrality')
barplot(scn.rel[order(scn.rel,decreasing=TRUE)],
        names=names(scn.rel)[order(scn.rel,decreasing=TRUE)],
        las=2,ylab='Relative Total Abundance')


##Tree scale
x.split <- paste(x$tree,x$quadrat,sep='_')
env.split <- paste(env$Tree.ID,env$Quad.Loc)
x.split <- as.character(x$tree)
env.split <- as.character(env$Tree.ID)
x.t <- split(x,x.split)
x.t <- split(x,x.split)
x.t <- lapply(x.t,function(x) x[,-c(1:4,ncol(x))])
cs <- unlist(lapply(x.t,cscore))
prb <- tapply(env$Pct.Roughness,env.split,mean) #percent rough bark
tid <- unlist(sapply(names(cs),function(x) strsplit(x,split='_')[[1]][1]))
                                        #co-occurrence nets
conets.t <- lapply(x.t,co.net)
net.d <- matrix(0,nrow=length(conets.t),ncol=length(conets.t))
rownames(net.d) <- colnames(net.d) <- names(conets.t)
for (i in 1:nrow(net.d)){
  for (j in 1:ncol(net.d)){
    net.d[i,j] <- sum(abs(conets.t[[i]]-conets.t[[j]])^2)
  }
}
                                        #
plot(cs~prb)
abline(lm(cs~prb))
summary(lm(cs~prb))
                                        #SES values
library(vegan)
library(pbapply)
#x.sim <- pblapply(x.t,nullCom)
#x.cs <- pblapply(x.sim,function(x) lapply(x,cscore))
#x.cs <- pblapply(x.cs,unlist)
                                        #dput(x.cs,file='~/projects/dissertation/projects/lichen_coo/results/x_cs.txt')
test <- dget(file='~/projects/dissertation/projects/lichen_coo/results/x_cs.txt')
x.cs <- test
ses <- cs*0
for (i in 1:length(ses)){
  ses[i] <- (cs[i] - mean(x.cs[[i]])) / sd(x.cs[[i]])
}
ses[is.na(ses)] <- 0
                                        #
par(mfrow=c(4,4))
for (i in 1:length(x.cs)){
  if (cs[i]<min(x.cs[[i]])){
    plot(density(x.cs[[i]]),main=names(cs)[i],xlab='C-Score',xlim=c(cs[i],max(x.cs[[i]])))
    abline(v=cs[i],lty=1,col='red')
  }else if (cs[i]>max(x.cs[[i]])){
    plot(density(x.cs[[i]]),main=names(cs)[i],xlab='C-Score',,xlim=c(min(x.cs[[i]]),cs[i]))
    abline(v=cs[i],lty=1,col='red')
  }
  else{
    plot(density(x.cs[[i]]),main=names(cs)[i],xlab='C-Score')
    abline(v=cs[i],lty=1,col='red')
  }
}

wild.p <- cs*0
for (i in 1:length(wild.p)){
  wild.p[i] <- length(x.cs[[i]][x.cs[[i]]<=cs[i]])/length(x.cs[[i]])
}

wild.coa <- cbind(ses=ses,p=wild.p)
write.csv(wild.coa,file='~/projects/dissertation/projects/lichen_coo/results/tables/wild_coa_table.csv')

###SES roughness, age, abundance, richness
height <- as.character(sapply(rownames(com),function(x) strsplit(x,split='_')[[1]][2]))
com.4555 <- com[height=='n45.55',]
com.8090 <- com[height=='n80.90',]
com.t <- com.4555+com.8090
total.abundance <- apply(com.t,1,sum)
species.richness <- apply(sign(com.t),1,sum)
                                        #
plot((ses~prb),xlab='Mean Percent Roughness',ylab='SES',pch=19,font.lab=2)
abline(lm(ses~prb))
summary(lm(ses~prb))
                                        #
plot((ses~total.abundance),xlab='Total Abundance',ylab='SES',pch=19,font.lab=2)
abline(lm(ses~total.abundance))
summary(lm(ses~total.abundance))
                                        #
plot((ses~species.richness),xlab='Species Richness',ylab='SES',pch=19,font.lab=2)
abline(lm(ses~species.richness))
summary(lm(ses~species.richness))
                                        #Age
age <- read.csv('~/projects/dissertation/projects/lichen_coo/data/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv')
                                        #predict age
gnu19.dbh <- age$DBH.cm_01[34]
new <- data.frame(DBH.cm_01=seq(min(age$DBH.cm_01),max(age$DBH.cm_01),by=0.1))
age <- na.omit(age)
pairs(cbind(age=age$AgeFinal.U,bole.circ=age$Bole.Circ.cm,DBH=age$DBH.cm_01,height=age$Height.m),lower.panel=panel.cor)
pred.age <- predict(lm(AgeFinal.U~DBH.cm_01,data=age),new)
plot(pred.age~new[,1])
gnu19.age <- as.numeric(pred.age[new[,1]==gnu19.dbh])
                                        #
age <- read.csv('~/projects/dissertation/projects/lichen_coo/data/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv')
age <- data.frame(tree.id=age[,1],age.final=age$AgeFinal.U)
age[,1] <- tolower(age[,1])
age[,1] <- sub('_','\\.',age[,1])
age[,1] <- sub('-','\\.',age[,1])
age[,1] <- sub('\\?','',age[,1])
age[,1] <- sub('\\.0','\\.',age[,1])
age[age[,1]=='gnu.85.1ftaway',1] <- 'gnu.85'
all(names(ses)%in%age[,1])
                                        #
tree.age <- numeric(length(ses))
tree.age <- age[match(names(ses),age[,1]),2]
tree.age[is.na(tree.age)] <- gnu19.age
                                        #
par(mfrow=c(1,2))
plot(prb~tree.age,xlab='Tree Age (years)',ylab='Percent Bark Roughness',font.lab=2,pch=19)
abline(lm(prb~tree.age))
summary(lm(prb~tree.age))
plot(ses~tree.age)
abline(lm(ses~tree.age))
summary(lm(ses~tree.age))
                                        #
library(sem)
###tree.age -> prb -> ses
model.sem <- specifyModel(file='~/projects/dissertation/projects/lichen_coo/src/sem_model.txt')
   sem.data <- cbind(tree.age,prb,ses)
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

###Microsite

env <- read.csv('~/projects/dissertation/projects/lichen_coo/data/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv')
env[,1] <- sub('\\?','',sub('\\.0','\\.',sub('\\_','\\.',sub('\\-','\\.',tolower(as.character(env[,1]))))))
env[env[,1]=='ll.6_(znu.29)',1] <- 'll.6'
env[env[,1]=='gnu.85.1ftaway',1] <- 'gnu.85'
env <- env[match(names(ses),env[,1]),]
range(na.omit(env$Height.m))

