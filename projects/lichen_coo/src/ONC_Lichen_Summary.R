###################################################
### chunk number 1: 
###################################################
#line 173 "ONC_Lichen_Summary.Rnw"
  ###Turn off the PerMANOVA's
  ###Make sure to turn this on when starting a new session in R.
  ###NOTE! You turned off the PerMANOVA for genotype alone
  adonis.on <- FALSE


###################################################
### chunk number 2: 
###################################################
#line 181 "ONC_Lichen_Summary.Rnw"
##Package dependencies
library(xtable)
library(ecodist)
library(vegan)
library(gdata)
library(gplots)
library(sna)
library(lme4)
library(RLRsim)
library(fossil)
source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')
source('/Users/Aeolus/Documents/R_Docs/Scripts/Functions/pairs.R')

##New functions
##Graph size (i.e. number of nodes)
netSize <- function(x){
  n <- apply(x,2,function(x) sum(abs(x)))
  n <- length(n[n!=0])
  return(n)
}

##Coerce tree notation
as.mytree <- function(x){

x <- as.character(x) #change from factor to character
x <- gsub('-','.',x) #change - to .

                                        #remove leading zeros in the tree number
for (i in (1:length(x))){
  y <- unlist(strsplit(x[i],split='\\.'))
  z <- unlist(strsplit(y[2],split=''))
  if (z[1] == '0'){
    x[i] <- paste(y[1],z[2],sep='.')
  }
}

return(x)
}




###################################################
### chunk number 3: 
###################################################
#line 226 "ONC_Lichen_Summary.Rnw"
##Load data
LCO <- read.csv('/Users/Aeolus/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')

##Remove genotypes
LCO <- LCO[LCO$Geno != 'T6',]
#LCO <- LCO[LCO$Geno != 'RL6',]
#LCO <- LCO[LCO$Geno != 'H10',]
#LCO <- LCO[LCO$Geno != '999',]
#LCO <- LCO[LCO$Geno != '1007',]

##Model networks for each quadrat
attach(LCO)
names(LCO)
tree.info <- paste(Tree,Geno,Quadrat,Year,sep='_') #integrate tree information
detach(LCO)

##Separate species observations
lco.com <- LCO[,7:15]



###################################################
### chunk number 4: 
###################################################
#line 250 "ONC_Lichen_Summary.Rnw"
##model networks
lco.nets <- list() #initiate tree network object

for (i in (1:length(unique(tree.info)))){
 x <- kendall.pairs(lco.com[tree.info == unique(tree.info)[i],],p.adj=FALSE)
 x[is.na(x)] <- 0#make NA values zero
 lco.nets[[i]] <- x
}

names(lco.nets) <- unique(tree.info) #name networks



###################################################
### chunk number 5: 
###################################################
#line 266 "ONC_Lichen_Summary.Rnw"
##Separate network factor information (tree, genotype, quadrat, year of sampling)
net.info <- unique(tree.info)
info.info <- rep(paste('tree','geno','quad','year',sep='_'),length(net.info))
tree.factor <- unlist(strsplit(net.info,split='_')) #tree-level information for each network
info.factor <- unlist(strsplit(info.info,split='_')) #sep information for tree information

net.tree <- tree.factor[info.factor == 'tree'] #tree number information for each network
net.geno <- tree.factor[info.factor == 'geno'] #genotype information for each network
net.quad <- tree.factor[info.factor == 'quad'] #quadrat information for each network
net.year <- tree.factor[info.factor == 'year'] #year of sampling for each network



###################################################
### chunk number 6: 
###################################################
#line 282 "ONC_Lichen_Summary.Rnw"
##Import roughness data
rough <- read.csv('data/ONC_raw_roughness.csv')
rough[is.na(rough)] <- 0 #replace NA values with zeros
rough[rough == ' '] <- 0 #replace empty spots with zeros
rough <- rough[rough[,1] != '',] #remove empty/place holder rows
names(rough)[1] <- 'Tree' #re-name the Tree column
names(rough)[2] <- 'Geno' #re-name the Genotype column
names(rough)[3] <- 'Quadrat' #re-name Quadrat column
                                        #coerce the Tree names
rough$Tree <- as.character(rough$Tree)
rough$Tree <- as.mytree(rough$Tree)
                                        #re-name quadrat locations
rough$Quadrat <- as.character(rough$Quadrat)
rough$Quadrat[rough$Quadrat == 'North 45-55'] <- 'n45.55' #north 45 - 55
rough$Quadrat[rough$Quadrat == 'North 80-90'] <- 'n80.90' #north 80 - 90
rough$Quadrat[rough$Quadrat == 'South 45-55'] <- 's45.55' #south 45 - 55
rough$Quadrat[rough$Quadrat == 'South 80-90'] <- 's80.90' #south 80 - 90

                                        #for each quadrat on each tree, separate the roughness information
names(rough)
attach(rough)
net.rough <- numeric()
rough.TQ <- paste(Tree,Quadrat) #create a vector for coalescing the roughness data
net.TQ <- paste(net.tree,net.quad)

for (i in (1:length(net.TQ))){
  net.rough[i] <- Roughness[rough.TQ == net.TQ[i]]
}

detach(rough)



###################################################
### chunk number 7: 
###################################################
#line 318 "ONC_Lichen_Summary.Rnw"
##Import tree position data
gps <- read.csv('data/NONC_gps_data.csv')[,1:3]
gps$Location <- as.mytree(as.character(gps$Location))
                                        #add missing trees
                                        #check with Jamie to make sure these data are correct
gps <- rbind(gps,c('N1.5', -111.9995,  41.24699))
gps <- rbind(gps,c('N1.7', -111.9995, 41.24709))
gps <- rbind(gps,c('N2.31', -111.9995823, 41.248322))
gps <- rbind(gps,c('N1.24', -111.999534, 41.247919))

                                        #initiate the lat and log vectors for the network locations
net.lat <- numeric() #lat
net.lon <- numeric() #lon

for (i in (1:length(net.tree))){
  net.lat[i] <- gps$Latitude[gps$Location == net.tree[i]]
  net.lon[i] <- gps$Longitude[gps$Location == net.tree[i]]
}

                                        #coerce into numeric
net.lat <- as.numeric(net.lat)
net.lon <- as.numeric(net.lon)

##Calculate species richness
net.rich <- numeric()
for (i in (1:length(unique(tree.info)))){
  x <- lco.com[tree.info == unique(tree.info)[i],]
  net.rich[i] <- sum(apply(x,2,function(x) if (any(x) > 0){1}else{0}))
}



###################################################
### chunk number 8: 
###################################################
#line 352 "ONC_Lichen_Summary.Rnw"
##PerMANOVA for networks
net.dist <- as.dist(netDist(lco.nets)) #calculate the network distances


                                        #Is space an important co-variate?
                                        #Probably NOT.
geo.dist <- earth.dist(cbind(net.lon,net.lat))
#mantel(geo.dist~net.dist)



###################################################
### chunk number 9: 
###################################################
#line 367 "ONC_Lichen_Summary.Rnw"
                                        #What about analyzing the quadrats separately?
                                        #YES, when we remove T6.
net.dist45 <- netDist(lco.nets[net.quad == 'n45.55'])
net.geno45 <- net.geno[net.quad == 'n45.55']
net.rough45 <- net.rough[net.quad == 'n45.55']
net.lon45 <- net.lon[net.quad == 'n45.55']
net.lat45 <- net.lat[net.quad == 'n45.55']
#adonis(net.dist45~net.geno45)
#adonis(net.dist45~net.geno45*net.rough45)
#adonis(net.dist45~net.geno45*net.lon45,permutations=9999)
Genotype <- factor(net.geno45)
Longitude <- net.lon45
#if (adonis.on == TRUE){adonis.genoXlon <- adonis(net.dist45~Genotype*Longitude,permutations=9999)}



###################################################
### chunk number 10: 
###################################################
#line 386 "ONC_Lichen_Summary.Rnw"
##Analyze structural statistics

##Only using the lower quadrats
lco.nets45 <- lco.nets[net.quad == 'n45.55']


##Degree
net.L45 <- unlist(lapply(lco.nets45,function(x) length(x[x!=0])/2)) #undirected, so only count one link per pair
#plot(net.L45~factor(net.geno45))
#summary(aov(net.L45~net.geno45*net.lon45))

##Size (i.e. number of nodes)
net.S45 <- unlist(lapply(lco.nets45,netSize))
#plot(net.S45~factor(net.geno45))
#plot(net.L45~net.S45)
#summary(aov(net.S45~net.geno45*net.lon45))

##Centralization
net.C45 <- unlist(lapply(lco.nets45,function(x) centralization(x,degree,mode='undirected')))
#hist(net.C45)
#hist(sqrt(asin(net.C45)))
#plot(net.C45~factor(net.geno45))
#summary(aov(sqrt(asin(net.C45))~net.geno45*net.lon45))


##Correlations with richness
net.rich45 <- net.rich[net.quad == 'n45.55']
Genotype <- factor(net.geno45)
Longitude <- net.lon45
aov.rich.genoXlon <- summary(aov(net.rich45~Genotype*Longitude))

#summary(aov(net.rich45~net.geno45*net.rough45))
#summary(aov(net.rich45~net.geno45*net.rough45))
#summary(aov(net.rich45~net.lat45))
#pairs(cbind(net.rich45,net.S45,net.L45,net.C45),pch=19)

##Whole model:
Richness <- net.rich45
Size <- net.S45
Degree <- net.L45
Centralization <- net.C45

if (adonis.on == TRUE){adonis.genoXrich <- adonis(net.dist45~Genotype*Richness,permutations=9999)}
aov.S45.genoXrich <- summary(aov(net.S45~Genotype*Richness))
aov.L45.genoXrich <- summary(aov(sqrt(net.L45)~Genotype*Richness))
aov.C45.genoXrich <- summary(aov(sqrt(asin(net.C45))~Genotype*Richness))
if (adonis.on == TRUE){adonis.netD.all <- adonis(net.dist45~Richness+Size+Degree+Centralization,permutations=9999)}



###################################################
### chunk number 11: gplots
###################################################
#line 440 "ONC_Lichen_Summary.Rnw"
##Graph plots
genomu.nets45 <- tapply(lco.nets45,net.geno45,edge.mean)
length(genomu.nets45)
par(mfrow=c(5,3),mai=c(0.3,0.3,0.3,0.3))
for (i in 1:length(genomu.nets45)){
  v.cex <- v.col<- apply(lco.com[LCO$Quadrat == 'n45.55' & LCO$Geno == names(genomu.nets45)[i],],2,sum)
  v.col[v.col > 0] <- 'black'
  v.col[v.col == 0] <- 'white'
  v.cex <- log(v.cex+1)
  v.cex <- v.cex + 1
  e.lwd <- (abs(genomu.nets45[[i]])*100)^9
  e.lwd[e.lwd > 25 & e.lwd < 50] <- 15
  e.lwd[e.lwd > 50] <- 25
  gplot(abs(genomu.nets45[[i]]),mode='circle',gmode='graph',vertex.sides=50,vertex.col=v.col,displaylabels=TRUE,vertex.border='black',vertex.cex=v.cex,label.cex=1.25,edge.lwd=e.lwd)
  title(main=names(genomu.nets45)[i],cex=5)
}




###################################################
### chunk number 12: ordinet
###################################################
#line 463 "ONC_Lichen_Summary.Rnw"
  ##Ordination of network statistics
pc.45 <- princomp(net.dist45)
env.45 <- cbind(net.S45,net.L45,net.rich45,net.C45,net.rough45,net.lon45)
ord.45 <- pc.45$scores[,1:2]
colnames(env.45) <- c('Size','Degree','Richness','Centralization','Roughness','Longitude')
vfit45 <- envfit(ord.45,env.45)

ord.mu45 <- cbind(tapply(ord.45[,1],net.geno45,mean),tapply(ord.45[,2],net.geno45,mean))
ord.se45 <- cbind(tapply(ord.45[,1],net.geno45,function(x) sd(x)/sqrt(length(x))),tapply(ord.45[,2],net.geno45,function(x) sd(x)/sqrt(length(x))))
ord.sd45 <- cbind(tapply(ord.45[,1],net.geno45,sd),tapply(ord.45[,2],net.geno45,sd))
ord.ciu45 <- ord.mu45 + ord.se45
ord.cil45 <- ord.mu45 - ord.se45

plot(ord.mu45,col=rainbow(max(as.numeric(factor(rownames(ord.mu45)))))[as.numeric(factor(rownames(ord.mu45)))],cex=2,xlim=c(min(ord.cil45[,1]),max(ord.ciu45)),ylim=c(min(ord.cil45[,2]),max(ord.ciu45[,2])),pch=as.numeric(factor(rownames(ord.mu45))),xlab='PC1',ylab='PC2')

for (i in 1:nrow(ord.mu45)){
  lines(c(ord.mu45[i,1],ord.mu45[i,1]),c(ord.cil45[i,2],ord.ciu45[i,2]),col=grey(0.25)[1])
  lines(c(ord.cil45[i,1],ord.ciu45[i,1]),c(ord.mu45[i,2],ord.mu45[i,2]),col=grey(0.25)[1])
}

plot(vfit45,col='black',lwd=5)
legend('bottomleft',legend=rownames(ord.mu45),pch=as.numeric(factor(rownames(ord.mu45))),col=rainbow(max(as.numeric(factor(rownames(ord.mu45)))))[as.numeric(factor(rownames(ord.mu45)))],bg='white',border='grey',cex=0.75)

#plot(pc.45,las=2)
#plot(pc.45$scores[,1:2],col=rainbow(max(as.numeric(factor(net.geno45))))[as.numeric(factor(net.geno45))],pch=19,cex=2)



###################################################
### chunk number 13: pairsplot
###################################################
#line 493 "ONC_Lichen_Summary.Rnw"
  p.plot <- cbind(net.rich45,net.S45,net.L45,net.C45)
colnames(p.plot) <- c('Richness','Size','Degree','Centralization')
pairs(p.plot,cex=0.85,upper.panel=panel.lm,lower.panel=panel.cor2)


###################################################
### chunk number 14: gplots
###################################################
#line 508 "ONC_Lichen_Summary.Rnw"
#line 440 "ONC_Lichen_Summary.Rnw#from line#508#"
##Graph plots
genomu.nets45 <- tapply(lco.nets45,net.geno45,edge.mean)
length(genomu.nets45)
par(mfrow=c(5,3),mai=c(0.3,0.3,0.3,0.3))
for (i in 1:length(genomu.nets45)){
  v.cex <- v.col<- apply(lco.com[LCO$Quadrat == 'n45.55' & LCO$Geno == names(genomu.nets45)[i],],2,sum)
  v.col[v.col > 0] <- 'black'
  v.col[v.col == 0] <- 'white'
  v.cex <- log(v.cex+1)
  v.cex <- v.cex + 1
  e.lwd <- (abs(genomu.nets45[[i]])*100)^9
  e.lwd[e.lwd > 25 & e.lwd < 50] <- 15
  e.lwd[e.lwd > 50] <- 25
  gplot(abs(genomu.nets45[[i]]),mode='circle',gmode='graph',vertex.sides=50,vertex.col=v.col,displaylabels=TRUE,vertex.border='black',vertex.cex=v.cex,label.cex=1.25,edge.lwd=e.lwd)
  title(main=names(genomu.nets45)[i],cex=5)
}


#line 509 "ONC_Lichen_Summary.Rnw"


###################################################
### chunk number 15: ordinet
###################################################
#line 528 "ONC_Lichen_Summary.Rnw"
#line 463 "ONC_Lichen_Summary.Rnw#from line#528#"
  ##Ordination of network statistics
pc.45 <- princomp(net.dist45)
env.45 <- cbind(net.S45,net.L45,net.rich45,net.C45,net.rough45,net.lon45)
ord.45 <- pc.45$scores[,1:2]
colnames(env.45) <- c('Size','Degree','Richness','Centralization','Roughness','Longitude')
vfit45 <- envfit(ord.45,env.45)

ord.mu45 <- cbind(tapply(ord.45[,1],net.geno45,mean),tapply(ord.45[,2],net.geno45,mean))
ord.se45 <- cbind(tapply(ord.45[,1],net.geno45,function(x) sd(x)/sqrt(length(x))),tapply(ord.45[,2],net.geno45,function(x) sd(x)/sqrt(length(x))))
ord.sd45 <- cbind(tapply(ord.45[,1],net.geno45,sd),tapply(ord.45[,2],net.geno45,sd))
ord.ciu45 <- ord.mu45 + ord.se45
ord.cil45 <- ord.mu45 - ord.se45

plot(ord.mu45,col=rainbow(max(as.numeric(factor(rownames(ord.mu45)))))[as.numeric(factor(rownames(ord.mu45)))],cex=2,xlim=c(min(ord.cil45[,1]),max(ord.ciu45)),ylim=c(min(ord.cil45[,2]),max(ord.ciu45[,2])),pch=as.numeric(factor(rownames(ord.mu45))),xlab='PC1',ylab='PC2')

for (i in 1:nrow(ord.mu45)){
  lines(c(ord.mu45[i,1],ord.mu45[i,1]),c(ord.cil45[i,2],ord.ciu45[i,2]),col=grey(0.25)[1])
  lines(c(ord.cil45[i,1],ord.ciu45[i,1]),c(ord.mu45[i,2],ord.mu45[i,2]),col=grey(0.25)[1])
}

plot(vfit45,col='black',lwd=5)
legend('bottomleft',legend=rownames(ord.mu45),pch=as.numeric(factor(rownames(ord.mu45))),col=rainbow(max(as.numeric(factor(rownames(ord.mu45)))))[as.numeric(factor(rownames(ord.mu45)))],bg='white',border='grey',cex=0.75)

#plot(pc.45,las=2)
#plot(pc.45$scores[,1:2],col=rainbow(max(as.numeric(factor(net.geno45))))[as.numeric(factor(net.geno45))],pch=19,cex=2)

#line 529 "ONC_Lichen_Summary.Rnw"


###################################################
### chunk number 16: pairsplot
###################################################
#line 547 "ONC_Lichen_Summary.Rnw"
#line 493 "ONC_Lichen_Summary.Rnw#from line#547#"
  p.plot <- cbind(net.rich45,net.S45,net.L45,net.C45)
colnames(p.plot) <- c('Richness','Size','Degree','Centralization')
pairs(p.plot,cex=0.85,upper.panel=panel.lm,lower.panel=panel.cor2)
#line 548 "ONC_Lichen_Summary.Rnw"


###################################################
### chunk number 17: 
###################################################
#line 565 "ONC_Lichen_Summary.Rnw"
  xtable(adonis.genoXrich$aov.tab,caption='ANOVA table for the PerMANOVA test of the effect of genotype and lichen species richness on the similarity of the lichen interaction networks.')


###################################################
### chunk number 18: 
###################################################
#line 574 "ONC_Lichen_Summary.Rnw"
  xtable(aov.S45.genoXrich)
  xtable(aov.L45.genoXrich)
  xtable(aov.C45.genoXrich,caption='ANOVA tables for the tests of the effects of genotype and richness on the network metrics (size, degree and centralization).')


###################################################
### chunk number 19: 
###################################################
#line 585 "ONC_Lichen_Summary.Rnw"
  xtable(adonis.netD.all$aov.tab,caption='ANOVA table for the PerMANOVA test of the effects of the network metrics (size, degree and centralization) on network similarity after removing the variance explained by richness.')


