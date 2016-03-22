###Generates a nested matrix from random draws from a binomial distribution

rnest <- function(dim=c(10,20),prnd=0.01){
  out <- rep(1,dim[1])
  for (i in 1:(dim[2]-1)){
    if (i <dim[1]){out <- c(out,sort(c(rep(0,i),rep(1,(dim[1]-i))),decreasing=TRUE))}else{out <- c(out,c(1,rep(0,(dim[1]-1))))}
  }
  out <- array(out,dim=dim)

  if (length(rnd)!=0){
    rnd <- sample(1:(dim[1]*dim[2]),round((length(1:(dim[1]*dim[2])))*prnd,1))
    out[rnd] <- abs(out[rnd]-1)
  }
  out
}
