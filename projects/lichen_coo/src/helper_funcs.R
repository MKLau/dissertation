
my.gplot <- function(scn,coord=''){
  if (coord==''){
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue')
}else{
  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue',coord=coord)
}
  
}
