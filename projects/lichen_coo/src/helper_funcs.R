
my.gplot <- function(scn){

  e.col <- sign(scn)
  e.col[e.col==1] <- 'grey'
  e.col[e.col==-1] <- 'red'
  gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue')
  
}
