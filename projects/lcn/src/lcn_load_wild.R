###LCN: loading the wild stand data
###MKLau
###21Mar2014

#library
library(ComGenR)
sapply(dir('~/projects/packages/ComGenR/R/',full.names=TRUE),source)
cs <- function(x){nestedchecker(x)[[1]][1]}
mm <- function(x){slot(computeModules(x),'likelihood')}
###
###
x <- read.csv('~/projects/dissertation/projects/lcn/data/lco_Apr2012.csv')
                                        #remove notes
x <- x[,colnames(x)!='NOTES.']
x <- x[,colnames(x)!='dead']
                                        #
x <- na.omit(x)
                                        #remove gnu.44 = FREMONT
x <- x[x$tree!='gnu.44',]
                                        #remove ll.6, weird tree with super smooth bark
x <- x[x$tree!='ll.6',]
x$tree <- factor(as.character(x$tree))
                                        #condense species
                                        #lecanora, there can be only one!
lec.sp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
                                        #no physcioids!
                                        #phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x,lec=lec.sp)
x <- x[,-c(6,8,10,18)]
x <- x[,colnames(x)!='physcioid']
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
colnames(x)[5:ncol(x)] <- c('Xg','Cs', 'Xm', 'fgb', 'Rs', 'Pm' ,'Pa', 'Pu','Ch','Ls')
x <- x[colnames(x)!='fgb']
x.q <- split(x,quads)
wild.com <- split(x,x$tree)
wild.com <- do.call(rbind,lapply(wild.com,function(x) apply(x[,-1:-4],2,sum)))
wild.q <- lapply(split(x,x$tree),function(x) x[,-1:-4])
                                        #data from lamit
env <- read.csv('~/projects/dissertation/projects/lcn/data/Uinta2012_all_data_from_Lamit.csv')
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
x.split <- paste(x$tree,x$quadrat,sep='_')
env.split <- paste(env$Tree.ID,env$Quad.Loc)
x.split <- as.character(x$tree)
env.split <- as.character(env$Tree.ID)
prb <- tapply(env$Pct.Roughness,env.split,mean) #percent rough bark
                                        #co-occurrence patterns
wco <- do.call(rbind,lapply(wild.q,function(x,t) apply(CoCo(x,type=t),2,sum),t='pos'))
wch <- do.call(rbind,lapply(wild.q,function(x,t) apply(CoCo(x,type=t),2,sum),t='neg'))
                                        #get ses values
                                        #"z" "means" "pval" "simulated" "method" "statistic" "alternative"
                                        #ws <- lapply(wild.q,function(x) oecosimu(x,cs,method='r1',burnin=100,thin=10,nsimul=5000))
                                        #wses <- unlist(lapply(ws,function(x) x$oecosimu[[1]]))
                                        #wsmu <- unlist(lapply(ws,function(x) x$oecosimu[[2]]))
                                        #wsp <- unlist(lapply(ws,function(x) x$oecosimu[[3]]))
                                        #wsim <- do.call(rbind,lapply(ws,function(x) x$oecosimu[[4]]))
                                        #rownames(wsim) <- names(wild.q)
                                        #ws <- cbind(wses,wsmu,wsp)
                                        #write.csv(ws,file='../data/wild_ses_21mar2014.csv')

###Rename data objects for simplicity
ws <- read.csv('~/projects/dissertation/projects/lcn/data/wild_ses_21mar2014.csv')
wc <- wild.com
wq <- wild.q
