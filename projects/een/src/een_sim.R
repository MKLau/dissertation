###Simulating community genetics patterns to study 
###the effects of foundation species genetic variation
###on species interactions via co-occurrence patterns.

library(ComGenR)
trees <- gpmTrees()
###com <- cgSim(tree.pheno=trees,reps=1,YY=30,GG=8)
com <- dget('data/cgsim_com.Rdata')
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

###Correlations with genetic variance
#get the genetic variance values for 

