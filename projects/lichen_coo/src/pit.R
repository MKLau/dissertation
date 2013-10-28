###LCO - Analysis of the Pit garden co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????


###Garden Analysis
garden.data <- read.csv('../data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
                                        #separate gardens
pit <- garden.data[garden=='pit',]
                                        #
pit.q <- split(pit,paste(pit$Tree,pit$Geno))
pit.tree <- unlist(sapply(names(pit.q),function(x) strsplit(x,split=' ')[[1]][1]))
pit.geno <- unlist(sapply(names(pit.q),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #remove 1012 for now
pit.q <- pit.q[pit.geno!='1012']
pit.tree <- pit.tree[pit.geno!='1012']
pit.geno <- pit.geno[pit.geno!='1012']
                                        #remove env data
pit.q <- lapply(pit.q,function(x) x[,-1:-6])
###Composition
com <- do.call(rbind,lapply(pit.q,function(x) apply(x,2,sum)))
com. <- cbind(com,ds=rep(1,nrow(com)))
com.rel <- apply(com.,2,function(x) x/max(x))
com.rel[is.na(com.rel)] <- 0
adonis(com.rel~pit.geno)
###SES analyses

