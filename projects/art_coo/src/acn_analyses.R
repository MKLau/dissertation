###Arthropod Co-occurrence Networks
###Modeling and analyses
###See notebook.Rnw for meta-data

binned.species <- TRUE

library(sna)
source('../../art_coo/src/helper_func.R')
source('../../lichen_coo/src/seenetR.R')

###Prelim data from ONC
## onc <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc.csv')
## onc.key <- read.csv('~/projects/dissertation/projects/art_coo/data/key.csv')
###Data from Pit

pit <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc_PIT_Lau.csv')
                                        #remove trailing 0
pit$tree <- sub('\\.0','\\.',pit$tree)
                                        #data verification of trees and genotypes
pg <- read.csv('~/projects/dissertation/docs/garden_information/PIT_garden_tree_information.csv')
pg <- pg[pg$Row!='DEAD',]
pg.tree <- as.character(pg$Row)
pg.tree <- sub('-','.',pg.tree)
pg.tree <- sub('NP','np',pg.tree)
                                        #genotype remove 1007
pit <- pit[pit$geno!='1007',]

###Data checks
                                        #na to zeros
pit[is.na(pit)] <- 0
                                        #merge categories
                                        #pemphigus mergers
pit$pb.upper <- pit$pb.upper + pit$pb.woody
pit$pb.pred <- pit$pb.pred + pit$pb.hole + pit$pb.woody.pred
pit <- pit[,colnames(pit)!='pb.woody'&colnames(pit)!='pb.woody.pred'&colnames(pit)!='pb.hole'&colnames(pit)!='mite']
                                        #remove species with less than 20 observations
pit.com <- pit[,-1:-6]
pit.com <- pit.com[,apply(pit.com,2,sum)>17]
                                        #separate live and senescing leaves
liv <- pit.com[pit[,1]=='live',]
sen <- pit.com[pit[,1]=='sen',]
                                        #make a list for each
liv.tl <- split(liv,pit$tree[pit[,1]=='live'])
sen.tl <- split(sen,pit$tree[pit[,1]=='sen'])
all(names(liv.tl)==names(sen.tl))
###Pemphigus patterns
                                        #pemphigus abundances
pb.liv <- unlist(lapply(liv.tl,function(x) sum(x$pb.upper)/length(x$pb.upper)))
pb.sen <- unlist(lapply(sen.tl,function(x) sum(x$pb.upper)/length(x$pb.upper)))
plot(pb.sen~pb.liv)
abline(lm(pb.sen~pb.liv))
hist((pb.sen-pb.liv))
t.test((pb.sen-pb.liv))
                                        #fungal
fung.liv <- unlist(lapply(liv.tl,function(x) sum(x$fungal)/length(x$fungal)))
fung.sen <- unlist(lapply(sen.tl,function(x) sum(x$fungal)/length(x$fungal)))
plot(fung.sen~fung.liv)
abline(lm(fung.sen~fung.liv))
hist((fung.sen-fung.liv))
t.test((fung.sen-fung.liv))

###Live leaves
###Stand level
par(mfrow=c(1,2))
liv.dnet <- dep.net(liv)
coord <- mgp2(liv.dnet,(apply(liv,2,sum)/max(apply(liv,2,sum)))+1,log.scale=FALSE)
sen.dnet <- dep.net(sen)
mgp2(sen.dnet,(apply(sen,2,sum)/max(apply(sen,2,sum)))+1,log.scale=FALSE,my.coord=coord)

###Tree level
                                        #separate trees
pit.l <- split(pit.com,paste(pit$tree,pit$geno))
tree <- as.character(unlist(sapply(names(pit.l),function(x) strsplit(x,split=' ')[[1]][1])))
geno <- as.character(unlist(sapply(names(pit.l),function(x) strsplit(x,split=' ')[[1]][2])))
                                        #co-occurrence
if (binned.species){
  pit.ses <- dget('../data/acn_tree_ses.Rdata')
  pit.pval <- dget('../data/acn_tree_pval.Rdata')
                                        #remove 1007
  g.ses <- as.character(unlist(sapply(names(pit.ses),function(x) strsplit(x,split=' ')[[1]][2])))
  pit.ses <- pit.ses[g.ses!='1007']
  pit.pval <- pit.pval[g.ses!='1007']
}else{
  pit.ses <- dget('../data/acn_tree_ses_nb.Rdata')
  pit.pval <- dget('../data/acn_tree_pval_nb.Rdata')
                                        #remove 1007
  g.ses <- as.character(unlist(sapply(names(pit.ses),function(x) strsplit(x,split=' ')[[1]][2])))
  pit.ses <- pit.ses[g.ses!='1007']
  pit.pval <- pit.pval[g.ses!='1007']
}
                                        #
pit.pval[pit.ses>0] <- 1 - pit.pval[pit.ses>0] 
pit.ses.zp <- pit.ses
pit.ses.zp[pit.pval>0.05] <- 0
hist(pit.ses)
hist(pit.ses.zp)
plot(pit.ses,pit.pval)
plot(pit.ses.zp,pit.pval)
summary(aov(pit.ses~geno))
par(mfrow=c(1,2))
plot(pit.ses~factor(geno))
                                        #plots
par(mfrow=c(1,2))
mgp2(pit.dnet,(apply(pit.com,2,sum)/max(apply(pit.com,2,sum)))+1,log.scale=FALSE)
plot(pit.ses~factor(geno),xlab='Genotype',ylab='SES')

                                        #genotype sensitivity
net.all <- pit.dnet
net.rnd <- dget(file='../data/acn_net_rnd.Rdata') #random removal stratified by genotype (5 trees)
net.grm <- dget(file='../data/acn_net_grm.Rdata') #genotype removal (number of trees equals genotype replicate)
                                        #
rnd.ncor <- numeric()
for (i in 1:length(net.rnd)){
  rnd.ncor[i] <- netCor(net.all,net.rnd[[i]])
}
                                        #
grm.ncor <- numeric()
for (i in 1:length(net.grm)){
  grm.ncor[[i]] <- netCor(net.all,net.grm[[i]])
}
hist(rnd.ncor)
abline(v=grm.ncor)
names(net.grm)[grm.ncor==min(grm.ncor)]
                                        #
rnd.ndist <- numeric()
for (i in 1:length(net.rnd)){
  rnd.ndist[i] <- sqrt(sum((net.all-net.rnd[[i]])^2))
}
                                        #
grm.ndist <- numeric()
for (i in 1:length(net.grm)){
  grm.ndist[[i]] <- sqrt(sum((net.all-net.grm[[i]])^2))
}
hist(rnd.ndist)
abline(v=grm.ndist)
names(net.grm)[grm.ndist==min(grm.ndist)]
grm.pval <- numeric()
for (i in 1:length(grm.ndist)){
  grm.pval[i] <- length(rnd.ndist[rnd.ndist>grm.ndist[i]])/length(rnd.ndist)
}

mark



##############################################################
############# OLD ANALYSES (before 15 Oct 2013) ##############
##############################################################

library(audio)
source('~/cor_nets/CorNets.R')
                                        #
                                        #data import
x <- read.csv('~/Dropbox/arth_cooc/arth_cooc_PIT_Lau.csv')
redo <- read.csv('~/Dropbox/arth_cooc/arth_cooc_PIT_redo.csv')
                                        #repalce NA with zeros
x[is.na(x)] <- 0
redo[is.na(redo)] <- 0
                                        #insert redos
x[x$tree=='np2.31',] <- redo[redo$tree=='np2.31',]
x[x$tree=='np4.36',] <- redo[redo$tree=='np4.36',]
all(x[x$tree=='np2.31',3:ncol(x)] == redo[redo$tree=='np2.31',3:ncol(x)])
all(colnames((redo[redo$tree=='np4.36',]))==colnames((x[x$tree=='np4.36',])))
                                        #remove live column for now
x <- x[,-1]
                                        #data checks
summary(x)
all(table(x$tree)==50)
table(x$geno)/50
                                        #find trees with pb
pb.trees <- unlist(lapply(split(x,x$tree),function(x) any(apply(x[,5:11],1,sum)!=0)))
names(pb.trees)[pb.trees]
                                        #find trees with pb.woody or pb.woody.pred
pbw.trees <- unlist(lapply(split(x,x$tree),function(x) any(apply(x[,8:9],1,sum)!=0)))
names(pbw.trees)[pbw.trees]
                                        #remove genotype 1007
x <- x[x$geno!='1007',]
x$geno <- factor(as.character(x$geno))
x$tree <- factor(as.character(x$tree))
table(x$geno)/50
table(x$tree)/50
                                        #combining modifiers
                                        #miners
miners <- apply(x[,16:18],1,sum)
x <- x[,-16:-18]
x <- cbind(x,miners)
                                        #combine pb and woody
pb <- apply(x[,c(5,8)],1,sum)
x[,5] <- pb
x <- x[,-8]
                                        #combine pb.pred and pb.woody.pred
pb.pred <- apply(x[,c(10,8)],1,sum)
x[,10] <- pb.pred
x <- x[,-8]
                                        #separte by tree
x.l <- split(x[,-1:-4],x$tree)
                                        #co-occurrence analysis by tree
car.l <- dget('./CAresults_trees')
## car.l <- lapply(x.l,CA.results) #car = co-occurrence analysis results
## play(sin(1:10000/10))
## car.bind <- do.call(rbind,lapply(car.l,function(x) x$cat))
## car.bind <- data.frame(car.bind)
## car.geno <- rownames(car.bind)
## for (i in 1:length(car.geno)){
##   car.geno[i] <- as.character(x$geno[x$tree==car.geno[i]][1])
## }

###differences in pb
pb <- apply(x[,5:9],1,sum)
pb <- tapply(pb,x$tree,sum)
summary(aov(pb~factor(car.geno)))
plot(pb~jitter(as.numeric(factor(car.geno)),0.5),xlab='Genotype',xaxt='n')
plot(pb~factor(car.geno),xlab='Genotype',xaxt='n',ylab=expression(paste(italic('P. betae'),' Abundance')))
axis(1,at=1:nlevels(factor(car.geno)),labels=levels(factor(car.geno)),las=2)

summary(fit <- aov(I(log(car.bind$SES+2))~factor(car.geno)))
shapiro.test(fit$residuals)
##                                         #
## par(mfrow=c(1,3))
## plot((I(log(car.bind$SES+2))~factor(car.geno)))
## plot((I((car.bind$SES+2))~factor(car.geno)))
## plot((I((car.bind$SES))~factor(car.geno)))
tukey.test <- TukeyHSD(fit)
tukey.test[[1]][tukey.test[[1]][,4]<=0.2,]

library(gplots)
mu.ses <- tapply(car.bind$SES,car.geno,mean)
se.ses <- tapply(car.bind$SES,car.geno,function(x) sd(x)/sqrt(length(x)))
barplot2(mu.ses,plot.ci=TRUE,ci.u=mu.ses+se.ses,ci.l=mu.ses-se.ses,las=2,ylab='SES',xlab='Genotype')
## plot(car.bind$SES~as.numeric(factor(car.geno)),xaxt='n',xlab='Genotype')
## axis(1,at=1:nlevels(factor(car.geno)),labels=levels(factor(car.geno)),las=2)
par(mfrow=c(2,2))
plot(car.bind$SES~pb,xlab='All pb types',ylab='SES')
abline(lm(car.bind$SES~pb))
summary(lm(car.bind$SES~pb))
                                        #
plot(car.bind$SES~tapply(x$pb.upper,x$tree,sum),xlab='pb.upper',ylab='SES')
abline(lm(car.bind$SES~tapply(x$pb.upper,x$tree,sum)))
summary(lm(car.bind$SES~tapply(x$pb.upper,x$tree,sum)))
                                        #
plot(car.bind$SES~tapply(x$pb.pred,x$tree,sum),xlab='pb.pred',ylab='SES')
abline(lm(car.bind$SES~tapply(x$pb.pred,x$tree,sum)))
summary(lm(car.bind$SES~tapply(x$pb.pred,x$tree,sum)))
                                        #
plot(car.bind$SES~tapply(x$pb.abort,x$tree,sum),xlab='pb.abort',ylab='SES')
abline(lm(car.bind$SES~tapply(x$pb.abort,x$tree,sum)))
summary(lm(car.bind$SES~tapply(x$pb.abort,x$tree,sum)))

  library(network)

                                        #co-occurrence networks
ca.nets <- lapply(t.com,function(x) araujoNet(x)$dp)
avg.nets <- list()
for (i in 1:length(unique(geno))){
  avg.nets[[i]] <- ca.nets[geno==unique(geno)[i]][[1]]
for (j in 2:length(ca.nets[geno==unique(geno)[i]])){
  x <- conform(avg.nets[[i]],ca.nets[geno==unique(geno)[i]][[j]])
  avg.nets[[i]] <- x[[1]] + x[[2]]
  }
  avg.nets[[i]] <- (avg.nets[[i]] / length(ca.nets[geno==unique(geno)[i]]))
}
                                        #convert to network objects
                                        #avg.nets <- lapply(avg.nets,as.network)

                                        #plot average networks
  par(mfrow=c(2,3),mar=c(5.1, 4.1, 4.1, 2.1)-2)
lapply(avg.nets,function(x) plot(as.network(abs(x)),edge.lwd=(abs(x)+2)^3,
                                 usearrows=FALSE,displaylabels=TRUE,vertex.cex=1.5,vertex.col='grey'))
