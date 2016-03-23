zeroCol <- function(x,n=10){
    for (i in 1:ncol(x)){
        if (sum(x[,i]) < n){x[,i] <- x[,i] * 0}else{}
    }
    x
}

coNets <- function(x){
    x <- sign(x)
    t(x) %*% x
}

rmZeros <- function(x,zero.diag=TRUE){
    if (zero.diag){diag(x) <- 0}
    x[apply(x,1,sum) != 0, apply(x,2,sum) != 0]
}

meanNet <- function(x,zero.diag=TRUE){
    mu <- x[[1]]
    for (i in 2:length(x)){
        mu <- mu + x[[i]]
    }
    mu <- mu / length(x)
    if (zero.diag){diag(mu) <- 0}
    return(mu)
}

varNet <- function(x,zero.diag=TRUE){
    mu <- meanNet(x,zero.diag=zero.diag)
    v <- (x[[1]] - mu)^2
    for (i in 2:length(x)){
        v <- v + (x[[i]] - mu)^2
    }
    v <- v / length(x)
    if (zero.diag){diag(v) <- 0}
    return(v)
}

distNet <- function(x,zero.diag=TRUE,dist.obj=TRUE){
    d <- matrix(0,nrow=length(x),ncol=length(x))
    for (i in 1:nrow(d)){
        for (j in 1:ncol(d)){
            xi <- x[[i]]
            xj <- x[[j]]
            if (zero.diag){diag(xi) <- diag(xj) <- 0}
            d[i,j] <- sqrt(sum((xi - xj)^2))
        }
    }
    if (dist.obj){d <- as.dist(d)}
    return(d)
}

ch.plot <- function(x='ordination matrix',g='groupings',cex=1,plot.legend=FALSE,loc='topleft',mu.pch=19){
    mu <- apply(x,2,function(x,g) tapply(x,g,mean),g=g) 
    se <- apply(x,2,function(x,g) tapply(x,g,function(x)
            sd(x)/sqrt(length(x))),g=g) 
    mu <- na.omit(mu) 
    se <- na.omit(se)
                                        #error bars
    cl.xu <- mu[,1] +  se[,1]
    cl.xl <- mu[,1] -  se[,1]
    cl.yu <- mu[,2] + se[,2]
    cl.yl <-  mu[,2] -  se[,2]
    if (plot.legend){
                                        #coloring
        mu.col <- rainbow(length(unique(g)))[as.numeric(unique(g))]
        plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
        for (i in 1:nrow(mu)){
            lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
            lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
        }
        legend(loc,legend=rownames(se),cex=cex*0.5,pch=mu.pch,col=mu.col,border='grey')
    }else{
                                        #coloring
        mu.col <- 'black'
        plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
        for (i in 1:nrow(mu)){
            lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
            lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
        }
    }
    mu
}
