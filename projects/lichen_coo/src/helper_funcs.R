
my.gplot <- function(scn,coord='',v.cex=''){
  if (v.cex==''){v.cex=rep(1,nrow(scn))}else{}
  if (coord==''){
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue',vertex.cex=v.cex)
}else{
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue',vertex.cex=v.cex
        ,coord=coord)
}
  
}
