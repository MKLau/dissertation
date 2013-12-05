###Arthropod Co-occurrence Networks
###Modeling and analyses
###See notebook.Rnw for meta-data

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
pb <- pit$pb.upper + pit$pb.lower
pit <- cbind(pit,pb=pb)
pit$pb.pred <- pit$pb.pred + pit$pb.hole + pit$pb.woody.pred
pit <- pit[,colnames(pit)!='pb.woody'&colnames(pit)!='pb.woody.pred'&colnames(pit)!='pb.hole'&colnames(pit)!='mite'&colnames(pit)!='pb.upper'&colnames(pit)!='pb.lower']
                                        #remove species with less than 17 observations
pit.com <- pit[,-1:-6]
pit.com <- pit.com[,apply(pit.com,2,sum)>17]
                                        #separate live and senescing leaves
liv <- pit.com[pit[,1]=='live',]
sen <- pit.com[pit[,1]=='sen',]
                                        #make a list for each
liv.tl <- split(liv,pit$tree[pit[,1]=='live'])
sen.tl <- split(sen,pit$tree[pit[,1]=='sen'])
all(names(liv.tl)==names(sen.tl))

###Live leaves
###Stand level
par(mfrow=c(1,2))
liv.dnet <- dep.net(liv)
coord <- mgp2(liv.dnet,(apply(liv,2,sum)/max(apply(liv,2,sum)))+1,log.scale=FALSE)
sen.dnet <- dep.net(sen)
mgp2(sen.dnet,(apply(sen,2,sum)/max(apply(sen,2,sum)))+1,log.scale=FALSE,my.coord=coord)

###Tree level
                                        #without fungus
ses.liv <- dget('../data/acn_ses_live_nf.Rdata')
ses.sen <- dget('../data/acn_ses_sen_nf.Rdata')
ses.liv[is.na(ses.liv)] <- 0
ses.sen[is.na(ses.sen)] <- 0
all(names(ses.liv)==names(ses.sen))
geno.ses <- unlist(sapply(names(ses.liv),function(x) strsplit(x,split=' ')[[1]][2]))
summary(aov(ses.liv~geno.ses))
summary(aov(ses.sen~geno.ses))
library(gplots)
par(mfrow=c(1,2))
mu <- tapply(ses.liv,geno.ses,mean)
se <- tapply(ses.liv,geno.ses,function(x) sd(x)/sqrt(length(x)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='SES',ylim=c(-3,3))
title(main='Living')
mu <- tapply(ses.sen,geno.ses,mean)
se <- tapply(ses.sen,geno.ses,function(x) sd(x)/sqrt(length(x)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='SES',ylim=c(-3,3))
title(main='Senescing')
                                        #geno effect of senescence
sen.eff <- ses.sen-ses.liv
summary(aov(sen.eff~geno.ses))
mu <- tapply(sen.eff,geno.ses,mean)
se <- tapply(sen.eff,geno.ses,function(x) sd(x)/sqrt(length(x)))
par(mfrow=c(1,1))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='(Sen - Liv) SES',ylim=c(-4,5))
title(main='Senescence Effect')
plot(TukeyHSD(aov(sen.eff~geno.ses)))
                                        #with fungus
ses.liv <- dget('../data/acn_ses_live.Rdata')
ses.sen <- dget('../data/acn_ses_sen.Rdata')
ses.liv[is.na(ses.liv)] <- 0
ses.sen[is.na(ses.sen)] <- 0
all(names(ses.liv)==names(ses.sen))
geno.ses <- unlist(sapply(names(ses.liv),function(x) strsplit(x,split=' ')[[1]][2]))
summary(aov(ses.liv~geno.ses))
summary(aov(ses.sen~geno.ses))
library(gplots)
par(mfrow=c(1,2))
mu <- tapply(ses.liv,geno.ses,mean)
se <- tapply(ses.liv,geno.ses,function(x) sd(x)/sqrt(length(x)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='SES',ylim=c(-5,2))
title(main='Living')
mu <- tapply(ses.sen,geno.ses,mean)
se <- tapply(ses.sen,geno.ses,function(x) sd(x)/sqrt(length(x)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='SES',ylim=c(-5,2))
title(main='Senescing')
                                        #geno effect on senescence effect
sen.eff <- ses.sen-ses.liv
summary(aov(sen.eff~geno.ses))
mu <- tapply(sen.eff,geno.ses,mean)
se <- tapply(sen.eff,geno.ses,function(x) sd(x)/sqrt(length(x)))
par(mfrow=c(1,1))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='(Sen - Liv) SES',ylim=c(-4,5))
title(main='Senescence Effect')
plot(TukeyHSD(aov(sen.eff~geno.ses)))

###Pemphigus patterns
                                        #pemphigus abundances
pb.liv <- unlist(lapply(liv.tl,function(x) sum(x$pb)/length(x$pb)))
pb.sen <- unlist(lapply(sen.tl,function(x) sum(x$pb)/length(x$pb)))
plot(pb.sen~pb.liv)
abline(lm(pb.sen~pb.liv))
hist((pb.sen-pb.liv))
t.test((pb.sen-pb.liv))
summary(aov(pb.liv~geno.ses))
summary(aov(pb.sen~geno.ses))
mu <- c(liv=mean(pb.liv),sen=mean(pb.sen))
se <- c(liv=sd(pb.liv)/sqrt(length(pb.liv)),sen=sd(pb.sen)/sqrt(length(pb.sen)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se)
pb.dif <- (pb.sen-pb.liv)
summary(aov(pb.dif~geno.ses))
mu <- tapply(pb.dif,geno.ses,mean)
se <- tapply(pb.dif,geno.ses,function(x) sd(x)/sqrt(length(x)))
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se)
                                        #
pbc.liv <- c(table(liv$pb[liv$pb!=0]),0)
names(pbc.liv)[4] <- "4"
pbc.sen <- table(sen$pb[sen$pb!=0])
pbc <- rbind(pbc.liv,pbc.sen)
fisher.test(pbc)
                                        #test of the number of twos
n2s <- tapply(pit$pb,paste(pit$tree,pit$leaf.type),function(x) length(x[x==2]))
n2s.leaf <- unlist(sapply(names(n2s),function(x) strsplit(x,split=' ')[[1]][2]))
summary(aov(n2s~n2s.leaf))
                                        #fungal
fung.liv <- unlist(lapply(liv.tl,function(x) sum(x$fungal)/length(x$fungal)))
fung.sen <- unlist(lapply(sen.tl,function(x) sum(x$fungal)/length(x$fungal)))
plot(fung.sen~fung.liv,pch=19)
abline(lm(fung.sen~fung.liv))
hist((fung.sen-fung.liv))
t.test((fung.sen-fung.liv))
summary(aov(fung.liv~geno.ses))
summary(aov(fung.sen~geno.ses))
summary(aov(I(fung.liv-fung.sen)~geno.ses))
