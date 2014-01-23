##Generate the stand level ses values
source('./seenetR.R')
garden.data <- read.csv('../data/LCO_data_ONC_PIT.csv')
                                        #remove genotype RL6 and N1.31
garden.data <- garden.data[garden.data$Geno!='RL6',]
garden.data <- garden.data[garden.data$Tree!='N1.31',]
                                        #separate onc
garden.data[,1] <- as.character(garden.data[,1])
g1 <- substr(garden.data[,1],2,2)
g1[g1!='P'] <- 'onc'
                                        #separate onc
onc <- garden.data[g1=='onc',]
					#stand
stand.null <- nullCom(onc[,7:15])
stand.null <- lapply(stand.null,cscore)
dput(stand.null,file='../data/onc_stand_null.Rdata')
