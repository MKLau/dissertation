### Re-making Fig2 for the final version of the manuscript
### 17Feb2016
### MKLau

setwd('~/Projects/pb_removal_nets/src/')
source('~/Projects/pb_removal_nets/src/pbr_load_data.R')
library(bipartite)

art <- pbr.gc08$c
x <- art;for (i in 1:20){print(c(i,ncol(x[,apply(art,2,sum) >= i])))}
art <- art[,apply(art,2,sum) >= 1]

plotweb(art)
