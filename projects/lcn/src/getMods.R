### getMods : takes a module object as input and returns the modules

getMods <- function(x,split.mods=TRUE){
  if (split.mods){
    com <- slot(x,'originalWeb')
    mod.s <- slot(x,'modules')[-1,-1:-2]
    mod.s <- apply(mod.s,2,function(x) (1:length(x))[x!=0])
    row.mods <- mod.s[1:nrow(com)]
    col.mods <- mod.s[(nrow(com)+1):(nrow(com)+ncol(com))]
    names(row.mods) <- rownames(com)
    names(col.mods) <- colnames(com)
    return(list(row.mods=row.mods,col.mods=col.mods))
  }else{
    mod.s <- slot(x,'modules')[-1,-1:-2]
    com <- slot(x,'originalWeb')
    names <- c(rownames(com),colnames(com))
    mod.l <- list()
    for (i in 1:nrow(mod.s)){mod.l[[i]] <- names[mod.s[i,mod.s[i,]!=0]]}
    return(mod.l)
  }
}
