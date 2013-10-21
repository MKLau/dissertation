##Generate the stand level ses values
source('./seenetR.R')

garden.data <- read.csv('../data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
                                        #separate onc
onc <- garden.data[garden=='onc',]
					#
stand.null <- nullCom(onc[,7:15])
stand.null <- lapply(stand.null,cscore)
dput(stand.null,file='../data/onc_stand_null.Rdata')
