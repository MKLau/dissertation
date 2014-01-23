###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????

rm(list=ls())
library(ComGenR)

###Garden Analysis
garden.data <- read.csv('~/projects/dissertation/projects/lcn/data/LCO_data_ONC_PIT.csv')
                                        #remove genotype RL6 and N1.31
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Tree!='N1.31',]
                                        #separate onc
garden.data[,1] <- as.character(garden.data[,1])
g1 <- substr(garden.data[,1],2,2)
g1[g1!='P'] <- 'onc'
onc <- garden.data[g1=='onc',]
					#tree overlap between years
unique(onc$Tree[onc$Year=='2010']) %in% unique(onc$Tree[onc$Year=='2011'])
unique(onc$Tree[onc$Year=='2011']) %in% unique(onc$Tree[onc$Year=='2010'])

###Microsat data from Nash
## library(polysat)
## library(xlsx)
## gen.d <- read.xlsx(file='../data/ONC_MSAT_datafromnash.xlsx',sheetIndex=2,row.names=TRUE)
## gen.d[is.na(gen.d)] <- 0
## gen.d <- as.dist((gen.d))
                                        #onc <- onc[onc$Year=='2011',]
###Merge species groups
#Phy <- apply(onc[,12:14],1,sum) #make phy out of pmel,pads,pund
#onc. <- onc[,-12:-14]
#onc <- data.frame(onc.,Phy)
if (all(table(onc[,1])==100)){}else{for (i in 1:1000){print('Warning: check input data!!!')}}
###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env);adonis(com~height);detach(env)

##Species accumulation curve
sac.com <- list()
for (i in 1:length(unique(env$tree))){
  sac.com[[i]] <- apply(com[env$tree==unique(env$tree)[i],],2,sum)
}
sac.com <- do.call(rbind,sac.com)
plot(specaccum(sac.com),ylim=c(0,10),ylab='Species Richness',xlab='Trees Sampled',font.lab=2)
legend('topright',legend=c('Garden','Wild'),lty=c(1,1),col=c(1,2))

###modeling
###Co-occurrences
##stand level
                                        #unipartite networks
lcn <- CoNetwork((onc[,7:ncol(onc)]))
mgp(lcn,(onc[,7:ncol(onc)]),displaylabels=TRUE,loc=FALSE)
                                        #degree and centrality vs abundance
deg <- apply(lcn,1,function(x) sum(sign(x)))
pca <- apply(onc[,7:ncol(onc)],2,function(x) sum(x)/length(x)) #percent abundance
cen <- sna::degree(lcn)
plot(cen~pca)
abline(lm(cen~pca))
summary(lm(cen~pca+I(pca^2)))
plot(deg~pca)
plot(log(deg+0.0001)~(pca))
summary(lm(log(deg+0.0001)~(pca)))
abline(lm(log(deg+0.0001)~(pca)))
                                        #bipartite network
                                        #genotype vector
tree.geno <- as.character(unlist(lapply(split(onc$Geno,onc$Tree),function(x) x[1])))
bpn.l <- split(onc[,7:ncol(onc)],onc$Tree)
bpn <- do.call(rbind,lapply(bpn.l,function(x) apply(x,2,sum)))
bpn.sort <- bpn[order(apply(bpn,1,function(x) sum(sign(x))),decreasing=TRUE),order(apply(bpn,2,function(x) sum(sign(x))),decreasing=TRUE)]
cgPlotweb(bpn,tree.geno,lab.cex=1.5,mean.geno=FALSE)
visweb(bpn.sort)
nest.r0 <- oecosimu(bpn,nestedtemp,method='r0',nsimul=1000,burnin=100)
nest.r1 <- oecosimu(bpn,nestedtemp,method='r1',nsimul=5000,burnin=100)
rich <- apply(bpn,1,function(x) sum(sign(x)))
shan <- diversity(bpn)
deg <- apply(bpn,2,function(x) sum(sign(x)))
toa <- apply(bpn,2,function(x) sum(x))
summary(aov(sqrt(shan)~factor(tree.geno)))
                                        #genotype network
bpn.geno <- mean.g(bpn,tree.geno)
nest.r1.geno <- oecosimu(bpn.geno,nestedtemp,method='r1',nsimul=1000,burnin=100)
cgPlotweb(bpn,tree.geno,lab.cex=2)
                                        #onc network
onc.cn <- co.net(onc[,7:ncol(onc)])
onc.dn <- dep.net(onc[,7:ncol(onc)])
onc.graph <- onc.dn[apply(onc.dn,1,sum)!=0,apply(onc.dn,2,sum)!=0]
par(mfrow=c(1,2))
my.gplot(onc.dn,v.cex=((apply(com,2,sum)/max(apply(com,2,sum)))+0.5))
onc.deg <- degree(onc.dn)
names(onc.deg) <- rownames(onc.cn)
barplot(onc.deg,ylab='Centrality')
                                        #co-occurrence
stand.null <- unlist(dget(file='../data/onc_stand_null.Rdata'))
stand.ses <- (cscore(onc[,7:ncol(onc)]) - mean(stand.null)) / sd(stand.null)
stand.ses.p <- length(stand.null[stand.null<=stand.ses])/length(stand.null)
c(stand.ses,stand.ses.p)
##tree level 
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
obs.cs <- unlist(lapply(onc.q,cscore))
                                        #load ses values
onc.ses <- dget('../data/onc_tree_ses.Rdata')
onc.p <- dget('../data/onc_tree_pval.Rdata')
ses.zero.p <- FALSE;ses.zero.sd2 <- FALSE
if (ses.zero.p){onc.ses[onc.p>0.05] <- 0}else{}
if (ses.zero.sd2){onc.ses[abs(onc.ses) < 2] <- 0}else{}
onc.tn <- lapply(onc.q,dep.net) #tree level networks
names(onc.ses) <- names(onc.q)[names(onc.q)!='']
onc.ses[is.na(onc.ses)] <- 0
onc.ses <- onc.ses[is.na(names(onc.ses))!=TRUE]
onc.cen <- unlist(lapply(onc.tn,function(x) centralization(x,FUN='degree')))
onc.deg <- unlist(lapply(onc.tn,function(x) length(x[x!=0])))
onc.netd <- netDist(onc.tn)
###Roughness in the Garden
rough <- read.csv('../data/ONC_raw_roughness.csv')
rough <- rough[as.character(rough[,1])!="",1:5]
                                        #isolate north quadrats
rough <- rough[sapply(rough[,3],function(x) substr(x,1,1)=='N'),]
                                        #average roughness
avg.rough <- tapply(rough[,5],rough[,1],mean)
r.tree <- names(avg.rough)
r.tree <- sub('-','\\.',r.tree)
r.tree <- sub('\\.0','\\.',r.tree)
names(avg.rough) <- r.tree
                                        #match roughness to to ses values
ses.tree <- as.character(sapply(names(onc.ses),function(x) unlist(strsplit(x,split=' '))[1]))
avg.rough <- avg.rough[match(ses.tree,r.tree)]
all(ses.tree==names(avg.rough))
###Genotype
genotype <- as.character(sapply(names(onc.q),function(x) unlist(strsplit(x,split=' '))[2]))

###Data match check
length(genotype)==length(avg.rough) & length(avg.rough)==length(onc.ses)
all(names(avg.rough)==as.character(unlist(sapply(names(onc.ses),function(x) strsplit(x,split=' ')[[1]][1]))))
all(names(onc.deg)==names(onc.ses))
###Analyses
                                        #community composition
onc.com <- do.call(rbind,lapply(onc.q,function(x) apply(x,2,sum)))
                                        #onc.com <- (onc.com/100)
onc.com <- cbind(onc.com,ds=rep(min(onc.com[onc.com!=0])/10,nrow(onc.com)))
adonis(onc.com~factor(genotype),permutations=10000)
adonis(onc.com~factor(genotype)+avg.rough,permutations=10000)
                                        #indicator species analysis
## library(labdsv)
## indval(onc.com[,-ncol(onc.com)],genotype)
## indval(onc.rel[,-ncol(onc.rel)],genotype)
## detach(package:labdsv)
                                        #what species like it rough
library(ecodist)
onc.nms <- nmds(vegdist(onc.com))
onc.ef <- envfit(nmds.min(onc.nms),env=avg.rough)
nms.col <- rainbow(max(as.numeric(factor(genotype))))[as.numeric(factor(genotype))]
plot(nmds.min(onc.nms),col=nms.col,pch=19,cex=2)
plot(onc.ef)
detach(package:ecodist)
                                        #networks and ses
summary(lm(onc.deg~onc.ses))
summary(lm(onc.cen~onc.ses))
summary(lm(onc.cen~onc.deg))
pairs(cbind(ses=onc.ses,degree=onc.deg,centrality=onc.cen))
                                        #roughness
summary(lm(log(abs(onc.ses)+0.1)~avg.rough))
summary(lm(onc.deg~avg.rough))
summary(lm(onc.cen~avg.rough))
adonis(onc.netd~avg.rough)
                                        #genotype
summary(aov(avg.rough~factor(genotype)))
plot(log(abs(onc.ses)+0.1)~factor(genotype))
summary(aov(log(abs(onc.ses)+0.1)~factor(genotype)))
summary(aov(log(abs(onc.ses)+0.1)~factor(genotype)))
summary(lm(log(abs(onc.ses)+0.1)~avg.rough))
shapiro.test(residuals(aov(log(abs(onc.ses)+0.1)~factor(genotype))))
hist(residuals(aov(log(abs(onc.ses)+0.1)~factor(genotype))))
ses.mu <- tapply(onc.ses,factor(genotype),mean)
ses.se <- tapply(onc.ses,factor(genotype),function(x) sd(x)/sqrt(length(x)))
barplot2(ses.mu,plot.ci=TRUE,ci.u=ses.mu+ses.se,ci.l=ses.mu-ses.se,las=2,ylab='SES',font.lab=2,cex.lab=1.25,font.axis=2)
                                        #ses predicts composition
library(ecodist)
mantel(vegdist(onc.com)~dist(onc.ses)+dist(avg.rough))
adonis(vegdist(onc.com)~onc.ses)
mantel(vegdist(onc.com)~dist(avg.rough)+dist(onc.ses))
mantel(vegdist(onc.com)~onc.netd)
mantel(onc.netd~dist(onc.ses))
                                        #nmds of composition
                                        #onc.nms <- nmds(vegdist(onc.com))
onc.nms <- dget(file='../results/figs/onc_nms.Rdata')
ch.plot(nmds.min(onc.nms),genotype,cex=2,plot.legend=FALSE,loc='topleft')
onc.ef <- envfit(nmds.min(onc.nms),env=data.frame(SES=onc.ses))
plot(onc.ef,labels='SES',col='black')
onc.ef
                                        #network correlations
                                        #onc.com <- com[,-ncol(com)]
wild.dn <- dget(file='~/projects/dissertation/projects/lcn/results/scn.Rdata')
wild.com <- dget(file='../data/wild_com.Rdata')
                                        #re-organize nodes to match
wild.dn <- wild.dn[c(1,2,6,5,3,4,7),c(1,2,6,5,3,4,7)]
wild.com <- wild.com[,c(1,2,6,5,3,4,7)]
rownames(wild.dn) <- colnames(wild.dn) <- colnames(onc.dn)
                                        #test for correlation in structure
par(mfrow=c(2,2))
coord <- mgp(wild.dn,wild.com)
mgp(onc.dn,onc.com,my.coord=coord)
plot(wild.dn[wild.dn!=0|onc.dn!=0]~onc.dn[wild.dn!=0|onc.dn!=0]
     ,xlab='Wild Network Edges',ylab='Garden Network Edges',pch=19,
     cex=1.5,font.lab=2,cex.lab=1.25,font.axis=2)
my.line(onc.dn[wild.dn!=0|onc.dn!=0],wild.dn[wild.dn!=0|onc.dn!=0],lwd=3)
cor.test(wild.dn[wild.dn!=0|onc.dn!=0],onc.dn[wild.dn!=0|onc.dn!=0])
delta.dn <- ((wild.dn-onc.dn)/onc.dn)
delta.dn[is.na(delta.dn)] <- 0
mgp2(abs(delta.dn),(abs(apply(wild.com,2,sum)-apply(onc.com,2,sum))/apply(onc.com,2,sum)),my.coord=coord,scalar=1)
                                        #pairs plot of species
pairs(onc.com[,-ncol(onc.com)])
                                        #genetic distance
onc.gavg <- matrix(0,nrow=length(unique(genotype)),ncol=ncol(onc.com))
rownames(onc.gavg) <- unique(genotype)
colnames(onc.gavg) <- colnames(onc.com)
for (i in 1:nrow(onc.gavg)){
  onc.gavg[i,] <- apply(onc.com[genotype==unique(genotype)[i],],2,mean)
  rownames(onc.gavg)[i] <- unique(genotype)[i]
}
                                        #match
gen.d <- as.matrix(gen.d)
gavg.d <- as.matrix(vegdist(onc.gavg))
gavg.d <- gavg.d[match(rownames(gen.d),rownames(gavg.d)),match(rownames(gen.d),rownames(gavg.d))]
gen.d <- as.dist(gen.d)
gavg.d <- as.dist(gavg.d)
                                        #mantel
mantel(gavg.d~gen.d)
                                        #ses
ses.mu <- tapply(onc.ses,genotype,mean)
ses.d <- as.matrix(dist(ses.mu))
ses.d <- as.dist(ses.d[match(rownames(as.matrix(gen.d)),rownames(ses.d)),match(rownames(as.matrix(gen.d)),rownames(ses.d))])
mantel(ses.d~gen.d)

