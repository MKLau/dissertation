zeroCol <- function(x,n=10){
    for (i in 1:ncol(x)){
        if (sum(x[,i]) < n){x[,i] <- x[,i] * 0}else{}
    }
    x
}
