###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????
library(sna)
source('./helper_funcs.R')

###Garden Analysis
garden.data <- read.csv('~/projects/dissertation/projects/lichen_coo/data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
                                        #separate gardens
onc <- garden.data[garden=='onc',]
pit <- garden.data[garden=='pit',]
                                        #onc network
par(mfrow=c(1,2))
onc.cn <- co.net(onc[,7:15])
onc.dn <- dep.net(onc[,7:15])
onc.graph <- onc.dn[apply(onc.dn,1,sum)!=0,apply(onc.dn,2,sum)!=0]
my.gplot(onc.graph)
                                        #
pit.cn <- co.net(pit[,7:15])
pit.dn <- dep.net(pit[,7:15])
pit.graph <- pit.dn[apply(pit.dn,1,sum)!=0,apply(pit.dn,2,sum)!=0]
my.gplot(pit.graph)

###Composition with height
library(vegan)
com <- split(onc[,7:ncol(onc)],paste(onc[,1],onc[,3],onc[,4]))
com <- do.call(rbind,lapply(com,function(x) apply(x,2,sum)))
com <- cbind(com,ds=rep(1,nrow(com)))
env <- data.frame(do.call(rbind,sapply(rownames(com),strsplit,split=' ')))
colnames(env) <- c('tree','year','height')
attach(env)

adonis(com~height)

detach(env)

###SES values
library(pbapply)
source('./seenetR.R')
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,3],onc[,4]))
onc.q <- lapply(onc.q,function(x) x[,7:15])
obs.cs <- unlist(lapply(onc.q,cscore))
onc.ses <- dget(file='~/projects/dissertation/projects/lichen_coo/data/onc_ses.Rdata')
## nonc.sim <- pblapply(onc.q,function(x) if (sum(sign(apply(x,2,sum)))>1){nullCom(x)}else{NA})
## onc.cs <- pblapply(onc.sim,function(x) if (any(is.na(x[[1]]))){NA}else{lapply(x,cscore)})
## onc.cs <- pblapply(onc.cs,unlist)
## onc.ses <- obs.cs*0
## for (i in 1:length(onc.ses)){
##   onc.ses[i] <- (obs.cs[i] - mean(onc.cs[[i]])) / sd(onc.cs[[i]])
## }


###Genotype and Roughness in the Garden
rough <- read.csv('~/projects/dissertation/projects/lichen_coo/data/ONC_raw_roughness.csv')
rough <- na.omit(rough[,1:5])
                                        #remove southern quadrats
rough[,3] <- as.character(rough[,3])
side <- sapply(rough[,3],function(x) unlist(strsplit(x,split=' '))[1])
rough <- rough[side=='North',]
                                        #reformat
rough[,1] <- as.character(rough[,1])
rough[,1] <- sub('\\-','\\.',rough[,1])
rough[rough[,3]=='North 45-55',3] <- 'n45.55'
rough[rough[,3]=='North 80-90',3] <- 'n80.90'
rough.tq <- paste(rough[,1],rough[,3])
rough.tq <- sub('\\.0','\\.',rough.tq)
onc.tq <- do.call(rbind,sapply(names(onc.ses),strsplit,split=' '))
onc.tq <- paste(onc.tq[,1],onc.tq[,3])
                                        #
onc.rough <- numeric(length(onc.tq))
for (i in 1:length(onc.tq)){
  if (onc.tq[i]%in%rough.tq){
    onc.rough[i] <- rough[rough.tq==onc.tq[i],5]
  }else{
    onc.rough[i] <- NA
  }
}

ses.data <- cbind(ses=onc.ses,rough=onc.rough)
plot(ses~rough,data=ses.data)
abline(lm(ses~rough,data=data.frame(ses.data)))
summary(lm(ses~rough,data=data.frame(ses.data)))

###Correlation of wild and onc network structure

###NEED TO RECONCILE TAXONOMIC DIFFERENCES
wild.net <- net.graph
onc.net <- onc.dn
colnames(wild.net)
colnames(onc.net)
