###LCN: loading the onc data
###MKLau
###21Mar2014

source('~/projects/packages/cooc/src/cooc.R')
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
                                        #
if (all(table(onc[,1])==100)){print('Good to go!')}else{for (i in 1:1000){print('Warning: check input data!!!')}}
                                        #separate trees
onc.q <- split(onc,paste(onc[,1],onc[,2]))
onc.q <- lapply(onc.q,function(x) x[,7:ncol(x)])
                                        #get genotype
onc.geno <- unlist(sapply(names(onc.q),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #Roughness in the Garden
rough <- read.csv('../../lcn/data/ONC_raw_roughness.csv')
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
get.ses <- function(x){print('...');if (all(x==0)){rep(NA,5)}else{oecosimu(x,cs,method='r1',burnin=100,thin=10,nsimul=5000)}}
onc.ses <- lapply(onc.q,get.ses)
os <- do.call(rbind,lapply(onc.ses,function(x) x[[2]][1:3]))
save(os,file='../data/lcn_onc_ses.rda')
