###LCO - Analysis of the wild stand co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????

###Tree data
x <- read.csv('~/dissertation/projects/lichen_coo/data/lco_Apr2012.csv')
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
env <- read.csv('~/dissertation/projects/lichen_coo/data/Uinta2012_all_data_from_Lamit.csv')
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
com.null <- nullCom(com)
com.obs <- cscore(com)
com.cs <- unlist(lapply(com.null,cscore))

##stand.null <- nullCom(x[,-c(1:4,ncol(x))])
stand.cs <- 
scn <- co.net(x[,-c(1:4,ncol(x))]) #stand-level co-occurrence network
gplot(abs(scn),displaylabels=TRUE,vertex.col='lightblue',gmode='graph',pad=1)



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
                      #
plot(cs~prb)
abline(lm(cs~prb))
summary(lm(cs~prb))
                                        #SES values
library(vegan)
library(pbapply)
x.sim <- pblapply(x.t,nullCom)
x.cs <- pblapply(x.sim,function(x) lapply(x,cscore))
x.cs <- pblapply(x.cs,unlist)
                                        #dput(x.cs,file='~/projects/dissertation/projects/lichen_coo/results/x_cs.txt')
test <- dget(file='~/dissertation/projects/lichen_coo/results/x_cs.txt')
x.cs <- test
hist(x.cs[[1]])
abline(v=cs[[1]])
ses <- cs*0
for (i in 1:length(ses)){
  ses[i] <- (cs[i] - mean(x.cs[[i]])) / sd(x.cs[[i]])
}
ses[is.na(ses)] <- 0
plot((ses~prb),xlab='Mean Percent Roughness',pch=19)
abline(lm(ses~prb))
summary(lm(ses~prb))
                                        #
par(mfrow=c(4,4))
for (i in 1:length(x.cs)){
  hist(x.cs[[i]],main=names(x.cs)[i],xlab='C-Score',)
  abline(v=cs[i],lty=1,col='red')
}
