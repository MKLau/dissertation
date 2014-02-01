###LCO - Analysis of the Pit garden co-occurrence data in Uintah, UT
###Taken out of the notebook.Rnw file chunk
###25 Sep 2013

###Meta

library(ComGenR)
library(lme4)
cgREML <- function(x,g){
  x.lmer <- lmer(x~(1|g))
  x.lm <- lm(x~1)
  chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
  p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  return(c(chi2=chi2,P.value=p.chi2))
}

###Garden Analysis
garden.data <- read.csv('../data/LCO_data_ONC_PIT.csv')
garden <- substr(garden.data[,1],2,2)
garden[garden=='P'] <- 'pit'
garden[garden!='pit'] <- 'onc'
                                        #separate gardens
pit <- garden.data[garden=='pit',]
                                        #
pit.q <- split(pit,paste(pit$Tree,pit$Geno))
pit.tree <- unlist(sapply(names(pit.q),function(x) strsplit(x,split=' ')[[1]][1]))
pit.geno <- unlist(sapply(names(pit.q),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #remove 1012 for now
pit.q <- pit.q[pit.geno!='1012']
pit.tree <- pit.tree[pit.geno!='1012']
pit.geno <- pit.geno[pit.geno!='1012']
                                        #remove env data
pit.q <- lapply(pit.q,function(x) x[,-1:-6])
###Composition
com <- do.call(rbind,lapply(pit.q,function(x) apply(x,2,sum)))
com. <- cbind(com,ds=rep(min(com[com!=0]),nrow(com)))
com.pc <- apply(com,2,function(x) x/length(x))
com.pc <- cbind(com.pc,ds=rep(min(com.pc[com.pc!=0],nrow(com.pc))))
com.rel <- apply(com,2,function(x) if (all(x==0)){x}else{x/max(x)})
com.rel <- cbind(com.rel,ds=rep(min(com.rel[com.rel!=0]),nrow(com)))
adonis(com.~pit.geno)
adonis(com.rel~pit.geno)
adonis(com.pc~pit.geno)

###Richness
R <- apply(sign(com),1,sum)
H <- diversity(com)
cgREML(R,pit.geno)
cgREML(H,pit.geno)

###SES analyses

###Nestedness
cgPlotweb(com,pit.geno)
pit.r00 <- oecosimu(com,nestfun=nestedtemp,method='r00',nits=1000)
pit.r0 <- oecosimu(com,nestfun=nestedtemp,method='r0',nits=1000)
pit.c0 <- oecosimu(com,nestfun=nestedtemp,method='c0',nits=1000)
pit.r1 <- oecosimu(com,nestfun=nestedtemp,method='r1',nits=1000)
pit.r00
pit.r0
pit.c0
pit.r1
