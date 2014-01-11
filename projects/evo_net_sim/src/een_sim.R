###Simulating community genetics patterns to study 
###the effects of foundation species genetic variation
###on species interactions via co-occurrence patterns.

library(ComGenR)
trees <- gpmTrees()
###com <- cgSim(tree.pheno=trees,reps=1,YY=30,GG=8)
com <- dget('../data/cgsim_com.Rdata')
geno <- factor(trees[,1])
test <- lapply(com,function(x,g) lapply(x,function(x,g) lapply(x,function(x,g) adonis(cbind(x,rep(1,nrow(x)))~g,permutations=1),g=g),g=g),g=geno)
test <- test[[1]]
r2 <- array(NA,dim=c(length(test),length(test[[1]])))
rownames(r2) <- paste('YY',1:length(test),sep='')
colnames(r2) <- paste('GG',1:length(test[[1]]),sep='')
for (i in 1:nrow(r2)){
  for (j in 1:ncol(r2)){
    r2[i,j] <- test[[i]][[j]]$aov.tab[1,5]
  }
}
net.l <- array(NA,dim=c(length(com[[1]]),length(com[[1]][[1]])))
rownames(net.l) <- paste('YY',1:length(com[[1]]),sep='')
colnames(net.l) <- paste('GG',1:length(com[[1]][[1]]),sep='')
for (i in 1:nrow(net.l)){
  for (j in 1:ncol(net.l)){
    com. <- com[[1]][[i]][[j]]
    com.[com.<=10] <- 0
    net <- CoNetwork(com.)
    net.l[i,j] <- length(net[net!=0])/2
  }
}
h2c <- array(NA,dim=c(length(com[[1]]),length(com[[1]][[1]])))
rownames(h2c) <- paste('YY',1:length(com[[1]]),sep='')
colnames(h2c) <- paste('GG',1:length(com[[1]][[1]]),sep='')
nms.dim <- h2c.sig <- r2.nms <- min.stress <- h2c
for (i in 1:nrow(h2c)){
  for (j in 1:ncol(h2c)){
    nms <- nmds(vegdist(com[[1]][[i]][[j]]),1,1,nits=5)
    min.nms <- nmds.min(nms)
    min.stress[i,j] <- nms$stress[nms$stress==min(nms$stress)][1]
    r2.nms[i,j] <- nms$r2[nms$stress==min(nms$stress)][1]
    h2c.sig[i,j] <- getH2C(min.nms,g=geno)[1]
    h2c[i,j] <- getH2C(min.nms,g=geno)[2]
    nms.dim[i,j] <- ncol(min.nms)
  }
}
print('Done!')
###Heatmaps
                                        #heritability significance
par(mfrow=c(1,2))
h2c. <- h2c
h2c.[h2c<=0.25] <- 0
image(t(h2c),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(h2c.),add=TRUE,nlevels=10,labcex=1,font.type=2)
title(main='H2C')
h2c.s <- h2c
h2c.s[h2c.sig<=0] <- 0
h2c.s. <- h2c.s
h2c.s.[h2c.sig<=0.25] <- 0
image(t(h2c.s),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(h2c.s.),add=TRUE,nlevels=10,labcex=1,font.type=2)
title(main='H2C.sig')
                                        #h2c vs r2
plot(as.vector(h2c)~as.vector(r2),xlab='R2P',ylab='H2C',font.lab=2)
abline(lm(as.vector(h2c)~as.vector(r2)))
h2r2.dif <- h2c-r2

par(mfrow=c(1,3))
plot(as.vector(h2r2.dif)~as.vector(min.stress),xlab='Stress',ylab='H2C - R2P',font.lab=2)
abline(lm(as.vector(h2r2.dif)~as.vector(min.stress)))
plot(as.vector(h2r2.dif)~as.vector(r2.nms),xlab='NMDS R2',ylab='H2C - R2P')
abline(lm(as.vector(h2r2.dif)~as.vector(r2.nms)))
plot(as.vector(r2.nms)~as.vector(min.stress),xlab='Stress',ylab='NMDS R2')
par(mfrow=c(1,1))
image(t(h2c),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(r2),add=TRUE,nlevels=10,labcex=1,font.type=2)
hist(I(t(h2c-r2)),font.lab=2,main='',xlab='H2C - R2permanova')
t.test(I(t(h2c-r2)))
                                        #
par(mfrow=c(1,3))
h2c. <- h2c
h2c.[h2c<=0.25] <- 0
image(t(h2c),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(h2c.),add=TRUE,nlevels=10,labcex=1,font.type=2)
title(main='H2C')
r2. <- r2
r2.[r2<=0.25] <- 0
image(t(r2),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(r2.),add=TRUE,nlevels=10,labcex=1,font.type=2)
title(main='PERMANOVA R2')
net.l. <- net.l/max(net.l)
image(t(net.l.),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(net.l.),add=TRUE,nlevels=20,labcex=1,font.type=2)
title(main='Network Degree')
                                        #
image(t(h2r2.dif),ylab='Ve',xlab='Vg',add=FALSE,font.lab=2)
contour(t(min.stress),add=TRUE,nlevels=20,labcex=1,font.type=2)
title(main='H2C - R2P')

###Quantiyfing Co-occurrence Vg

cg.test <- function(x='community matrix',g='grouping factor',fdr=TRUE,co.return=FALSE){
  ##Get co-occurrences
  x[x!=0] <- 1
  co.mat <- list()
  k <- 1
  for (i in 1:(ncol(x))){
    for (j in i:ncol(x)){
      if (i == ncol(x)|i==j){}else{
        co.mat[[k]] <- x[,i] + x[,j]
        names(co.mat)[k] <- paste(i,j,sep='_')
        k <- k + 1
      }
    }
  }
  co.mat <- do.call(cbind,co.mat)
  co.mat[co.mat!=2] <- 0
  co.mat[co.mat==2] <- 1
  if (co.return){return(co.mat)}else{
    co.test <- apply(co.mat,2,function(x,g) glm(x~g,family='binomial'),g=g)
    p.val <- unlist(lapply(co.test,function(x) unlist(anova.glm(x,test='Chisq'))[10]))
    names(p.val) <- names(co.test)
    if (fdr){pval <- p.adjust(p.val,method='fdr')}
    return(p.val)
  }
}

abund.shift <- function(x){
  shift.m <- list()
  k <- 1
  for (i in 1:(ncol(x))){
    for (j in i:ncol(x)){
      if (i == ncol(x)|i==j){}else{
        shift.m[[k]] <- abs(x[,i]-x[,j]) / apply(cbind(x[,i],x[,j]),1,max)
        names(shift.m)[k] <- paste(i,j,sep='_')
        k <- k + 1
      }
    }
  }
  shift.m <- do.call(cbind,shift.m)
  shift.m[is.na(shift.m)] <- 0
  return(shift.m)
}

my.com <- com[[1]][[1]][[7]]
my.com[my.com<=10] <- 0
co.mat <- cg.test(my.com,factor(geno),fdr=TRUE,co.return=TRUE)
a.mat <- abund.shift(my.com)
adonis(my.com~g)
adonis(co.mat~g)
adonis(a.mat~g)

am.test <- apply(a.mat,2,function(x,g) glm(x~g),g=factor(geno))
p.val <- unlist(lapply(am.test,function(x) unlist(anova.glm(x,test='F'))[12]))
p.val <- p.adjust(p.val)
names(p.val) <- names(am.test)
names(p.val)[p.val<=0.05]
am.net <- sapply(names(p.val)[p.val<=0.05],function(x) strsplit(x,split='_')[[1]])
am.am <- array(0,dim=c(length(unique(as.vector(am.net))),length(unique(as.vector(am.net)))))
rownames(am.am) <- colnames(am.am) <- unique(as.vector(am.net))
for (i in 1:ncol(am.net)){
  am.am[rownames(am.am)==am.net[1,i],colnames(am.am)==am.net[2,i]] <- 1
  am.am[rownames(am.am)==am.net[2,i],colnames(am.am)==am.net[1,i]] <- 1
}
gplot(am.am)

am.test <- apply(a.mat,2,function(x,g) glm(x~g),g=factor(geno))
