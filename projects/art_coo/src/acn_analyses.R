###Arthropod Co-occurrence Networks
###Modeling and analyses
###See notebook.Rnw for meta-data

library(sna)
source('../../lichen_coo/src/helper_funcs.R')
source('../../lichen_coo/src/seenetR.R')

###Prelim data from ONC
## onc <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc.csv')
## onc.key <- read.csv('~/projects/dissertation/projects/art_coo/data/key.csv')
###Data from Pit

pit <- read.csv('~/projects/dissertation/projects/art_coo/data/arth_cooc_PIT_Lau.csv')

###Data checks
                                        #na to zeros
pit[is.na(pit)] <- 0
                                        #merge categories
pit.com <- pit[,-c(1,2,3,4,5)]
pb <- apply(pit.com[,1:7],1,sum)
pit.com <- cbind(pb,pit.com[,-1:-7])
pit.com <- cbind(pit.com[,-c(2,3)],chew=apply(pit.com[,2:3],1,sum))
pit.com <- pit.com[,apply(pit.com,2,sum)>1]
pit.com$pb[pit.com$pb!=0] <- 1
pit.com$chew[pit.com$chew!=0] <- 1
pit.com <- as.matrix(pit.com)

###Stand level
pit.dnet <- dep.net(pit.com)
mgp2(pit.dnet,(apply(pit.com,2,sum)/max(apply(pit.com,2,sum)))+1,log.scale=FALSE)

###Tree level
                                        #separate trees
pit.l <- split(pit.com,paste(pit$tree,pit$geno))
geno <- as.character(unlist(sapply(names(pit.l),function(x) strsplit(x,split=' ')[[1]][2])))
                                        #co-occurrence




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
