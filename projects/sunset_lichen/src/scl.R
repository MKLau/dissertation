###Sunset Lichen Networks

key <- read.csv('../data/2012jul30/Key.csv')
x <- read.csv('/Users/Aeolus/projects/dissertation/projects/sunset_lichen/data/2012jul30/spp_env_combined.csv')
                                        #remove dead
x <- x[x$Live.Dead == 1,]
com <- x[,((1:ncol(x))[colnames(x) == 'Acacon']):((1:ncol(x))[colnames(x) == 'Xanele'])]
env <- x[,1:12]
                                        #remove N and S light
env <- env[,-10:-11]
                                        #fix colnames
colnames(env) <- sub('\\.\\.','',colnames(env))
colnames(env)
                                        #bipartite network
env. <- env[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE),]
bpn <- com[order(apply(com,1,function(x) sum(sign(x))),decreasing=TRUE),order(apply(com,2,function(x) sum(sign(x))),decreasing=TRUE)]
rownames(bpn) <- paste(env.$Tree.pairs,env.$Moth,sep='_')
plotweb(bpn,method='normal',text.rot=90,col.low=env.$Moth+1)
nestedness(bpn)
