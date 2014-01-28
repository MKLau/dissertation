###Getting SES values for een_sym

library(vegan)
library(bipartite)
###Sourcing ComGenR on hoth
oldwd <- getwd()
setwd('../../../../ComGenR/R/')
cgn.list <- (sapply(dir(),grepl,pattern='~')|sapply(dir(),grepl,pattern='\\#'))==FALSE
sapply(dir()[cgn.list],source)
setwd(oldwd)
###
x <- dget('../results/een_exp_sym/comsim_sym.rdata')
out <- lapply(x,function(x) lapply(x,function(x) cnm.test(x,nits=5000)[1]))
names(out) <- paste('sim',(1:length(x)),sep='')
write.csv(do.call(rbind,out),file='../results/een_exp_sym/een_ses.csv')
