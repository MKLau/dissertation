

###cross hair plots

ch.plot <- function(x,g,cex=1,buffer=0.1,plot.legend=TRUE,loc='topleft'){
  mu <- apply(x,2,function(x,g) tapply(x,g,mean),g=g)
  se <- apply(x,2,function(x,g) tapply(x,g,function(x) sd(x)/sqrt(length(x))),g=g)
  mu <- na.omit(mu)
  se <- na.omit(se)

                                        #coloring
  mu.col <- 'black'
  mu.pch <- 19

                                        #error bars
  cl.xu <- mu[,1] + se[,1]
  cl.xl <- mu[,1] - se[,1]
  cl.yu <- mu[,2] + se[,2]
  cl.yl <- mu[,2] - se[,2]
  plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
  for (i in 1:nrow(mu)){
    lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
    lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
  }
  if (plot.legend){legend(loc,legend=rownames(se),cex=cex*0.5,pch=mu.pch,col=mu.col,border='grey')}else{}
}

## Usage
## xtype <- as.character(x.[,colnames(x.)=="Genotype..FINAL."])
## xtype[xtype=='BC1-DxAA '] <- 'BC1-DxAA'
## xtype[xtype=='trihybride'] <- 'Trihybrid'
## ch.plot(nmds.min(nms.),xtype,cex=2,plot.legend=TRUE,loc='topright')


mgp <- function(scn,com,loc=TRUE,my.coord=''){
e.col <- sign(scn)
e.col[e.col==1] <- 'grey'
e.col[e.col==-1] <- 'red'
v.cex <- apply(com,2,sum)
v.cex <- log(v.cex,10)
v.cex <- v.cex * (1/min(v.cex))
v.cex <- v.cex/2
if (length(my.coord)==1){
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey')
}else{
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey',
                 coord=my.coord)
}
if (loc){return(coord)}else{}
}

mgp2 <- function(scn,v.cex,loc=TRUE,my.coord='',log.scale=TRUE,scalar=2){
e.col <- sign(scn)
e.col[e.col==1] <- 'grey'
e.col[e.col==-1] <- 'red'
if (log.scale){
  v.cex <- log(v.cex,10)
  v.cex <- v.cex * (1/min(v.cex))
  v.cex <- v.cex/scalar
  is.na(v.cex) <- 0
}else{}
if (length(my.coord)==1){
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey',)
}else{
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey',
                 coord=my.coord)
}
if (loc){return(coord)}else{}
}

my.line <- function(x,y,lwd=2){
  fit <- lm(y~x)
  lines(x,fitted.values(fit),lwd=lwd)
}

my.gplot <- function(scn,coord='',v.cex=''){
  if (v.cex==''){v.cex=rep(1,nrow(scn))}else{}
  if (coord==''){
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='darkgrey',vertex.cex=v.cex)
}else{
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='darkgrey',vertex.cex=v.cex
        ,coord=coord)
}
  
}

netDist <- function(dn.t){
  net.d <- matrix(0,nrow=length(dn.t),ncol=length(dn.t))
  rownames(net.d) <- colnames(net.d) <- names(dn.t)
  for (i in 1:nrow(net.d)){
    for (j in 1:ncol(net.d)){
      net.d[i,j] <- sum(abs(dn.t[[i]]-dn.t[[j]])^2)
    }
  }
  net.d <- as.dist(net.d)
  return(net.d)
}

netCor <- function(net1,net2){
  ncor <- cor(net1[net1!=0|net2!=0],net2[net1!=0|net2!=0])
}
