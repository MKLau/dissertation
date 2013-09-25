###LCO - Analysis of the ONC garden co-occurrence data in Ogden, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta
##?????


###Garden Analysis
garden.data <- read.csv('~/projects/dissertation/projects/onc_lichen/data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
                                        #separate gardens
onc <- garden.data[garden=='onc',]
pit <- garden.data[garden=='pit',]
                                        #onc network
onc.cn <- co.net(onc[,7:15])
onc.dn <- dep.net(onc[,7:15])
onc.graph <- onc.dn[apply(onc.dn,1,sum)!=0,apply(onc.dn,2,sum)!=0]
gplot(onc.graph,displaylabels=TRUE,edge.col='darkgrey',
      gmode='graph',label.cex=1.5,vertex.col='lightblue',
      edge.lwd=(onc.graph/max(onc.graph)+1)^2.5,uselen=FALSE,
      edge.len=0.025,usearrows=FALSE)
                                        #pit network
pit.cn <- co.net(pit[,7:15])
pit.dn <- dep.net(pit[,7:15])
pit.graph <- pit.dn[apply(pit.dn,1,sum)!=0,apply(pit.dn,2,sum)!=0]
gplot(pit.graph,displaylabels=TRUE,edge.col='darkgrey',
      gmode='graph',label.cex=1.5,vertex.col='lightblue',
      edge.lwd=(pit.graph/max(pit.graph)+1)^2.5,uselen=FALSE,
      edge.len=0.025,usearrows=FALSE)

###Roughness in the Garden
rough <- read.csv('data/ONC_raw_roughness.csv')
head(rough)

###Correlation of wild and onc network structure

###NEED TO RECONCILE TAXONOMIC DIFFERENCES
wild.net <- net.graph
onc.net <- onc.dn
colnames(wild.net)
colnames(onc.net)
