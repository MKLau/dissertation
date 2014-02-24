###Modularity analysis for the een networks
###M.K.Lau 24 Apr 2014

args <- commandArgs(trailingOnly = TRUE)

library(bipartite)
library(methods)
file.path <- as.character(args[1])

print(file.path)
x <- read.csv(file.path)
m <- sign(x)[,apply(x,2,sum)!=0]
print('Null Modules')
                                        #burning
null <- commsimulator(m,method='r1',thin=100)
for (i in 1:500){null <- commsimulator(null,method='r1',thin=100)}
                                        #
null <- commsimulator(null,method='r1',thin=100)
mod.null <- slot(computeModules(null),name='likelihood')
cat(mod.null,file='../results/bitmodout.txt',append=TRUE,sep=' ')


