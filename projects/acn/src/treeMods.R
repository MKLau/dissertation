
tree.mod <- list()
for (i in 45:length(tree.nets)){
    if (sum(sign(apply(tree.nets[[i]],2,sum))) < 3){
        tree.mod[[i]] <- 0
    }else{
        out <- computeModules(tree.nets[[i]])
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
