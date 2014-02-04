###LCO - Analysis of the wild stand co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###4 Feb 2014

###Meta
##Site = Uintah, UT
##Study area = 225 * 463 = 104,175 m2 = 0.104175 km2

library(pbapply)
library(vegan)
build.bpn <- function(x,alpha=0.05,p=0.001,adjust=FALSE){
  p.out <- apply(x,2,function(x) as.numeric(unlist(binom.test(sum(sign(x)),length(x),p=p))[3]))
  if (adjust){p.adjust(p.out,method='fdr')}
  x.out <- apply(sign(x),2,sum)/nrow(x)
  x.out[p.out>alpha] <- 0
  return(x.out)
}

                                        #Sourcing ComGenR on hoth
oldwd <- getwd()
setwd('~/projects/packages/ComGenR/R/')
cgn.list <- (sapply(dir(),grepl,pattern='~')|sapply(dir(),grepl,pattern='\\#'))==FALSE
sapply(dir()[cgn.list],source)
setwd(oldwd)

###lichen data
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
lec.sp <- apply(x[,c(6,8,10,18)],1,function(x) sign(any(x!=0)))
                                        #phy.spp <- apply(x[,c(13,14,15,16)],1,function(x) sign(any(x!=0)))
x <- cbind(x,lec=lec.sp)
x <- x[,-c(6,8,10,18)]
                                        #break into quadrat list (x.q)
quads <- paste(x$tree,x$quadrat)
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
                                        #environmental data from lamit
wenv <- read.csv('~/projects/dissertation/projects/lcn/data/uintah_trees/uintah_2012_lco_trees.csv')
wenv$tree <- sub('\\.0','\\.',wenv$tree) #remove leading zeros
wenv$tree[wenv$tree=='gnu.468'] <- "gnu.68.ramet"
                                        #match to lichen dataset
wenv <- wenv[wenv$tree%in%names(prb),]
wenv <- wenv[match(names(prb),wenv$tree),]
                                        #lat and long
lat <- wenv$lat..N.deg_decimal.
lat <- t(sapply(lat,function(x) as.numeric(unlist(strsplit(as.character(x),split='_')))))
lat <- lat[,1] + lat[,2]/60
lon <- wenv$lon
lon <- t(sapply(lon,function(x) as.numeric(unlist(strsplit(as.character(x),split='_')))))
lon <- lon[,1] + lon[,2]/60
wenv$lat <- lat
wenv$lon <- lon
                                        #calculate geographic distance
library(fossil)
geo.dist <- earth.dist(data.frame(lon,lat))
range(geo.dist * 1000)

#roughness and community composition
wild.com.rel <- apply(wild.com,2,function(x) if(all(x==0)){x}else{x/max(x)})
adonis(wild.com~prb)
adonis(wild.com.rel~prb)

#ses patterns
print('stand ses')
wild.stand.ses <- cnm.test(wild.com)
print(wild.stand.ses)
print('tree ses')
wild.tree.ses <- lapply(wild.q,cnm.test,nits=5000)
wild.ses <- do.call(rbind,wild.tree.ses)[,1]
wild.ses[is.na(wild.ses)] <- 0
write.csv(do.call(rbind,wild.tree.ses),file='../../lcn/data/wild_ses_tree.csv')

#nestedness
print('nestedness')
wild.bpn <- do.call(rbind,lapply(wild.q,build.bpn))
                                        #cgPlotweb(wild.bpn,rownames(wild.bpn))
wild.nest <- list()
wild.nest[[1]] <- oecosimu(wild.bpn,nestfun=nestedtemp,method='r00',nsimul=5000)
wild.nest[[2]] <- oecosimu(wild.bpn,nestfun=nestedtemp,method='r0',nsimul=5000)
wild.nest[[3]] <- oecosimu(wild.bpn,nestfun=nestedtemp,method='c0',nsimul=5000)
wild.nest[[4]] <- oecosimu(wild.bpn,nestfun=nestedtemp,method='r1',nsimul=5000)
wild.nest.out <- do.call(rbind,lapply(wild.nest,function(x) x$oecosimu[-4]))
write.csv(wild.nest.out,file='../../lcn/results/wild_nest.csv')
print('done.')
