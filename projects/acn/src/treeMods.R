
tree.mod <- list()
for (i in 1:length(tree.nets)){
    x <- tree.nets[[i]]
    diag(x) <- 0
    if (sum(sign(apply(x,2,sum))) < 3){
        tree.mod[[i]] <- 0
    }else{
        out <- computeModules(x)
        tree.mod[[i]] <- slot(out,'likelihood')
    }
}
tree.mod <- unlist(tree.mod)
dput(tree.mod,file='../data/tree.mod')

### live and senescing bipartite networks
liv.modules <- computeModules(liv.bpn)
sen.modules <- computeModules(sen.bpn)
dput(liv.modules,'../data/liv.modules')
dput(sen.modules,'../data/sen.modules')
