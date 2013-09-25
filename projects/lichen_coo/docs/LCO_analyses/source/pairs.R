## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * (abs(r)+0.2) )
#    text(0.5, 0.5, txt, cex = 2)
}

#Resize
panel.kendall <- function(x, y, digits=2, prefix="", cex.cor,method='kendall')
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,method=method)
    p.cor<-cor.test(x,y,method=method,alternative='two.sided',exact=FALSE)
    p.cor<-as.numeric(unlist(p.cor)[2])
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    cex.cor.r <-cex.cor * abs(r)
    if (p.cor<=0.05){text(0.5, 0.5, txt, cex = 1,col='blue')}else {text(0.5, 0.5, txt, cex = 0.75)}

}

#Bonferroni Correction
panel.bfcor <- function(x, y, digits=2, prefix="", cex.cor,method='spearman',alpha=0.1)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,method=method)
    p.cor<-cor.test(x,y,method=method,alternative='two.sided',exact=FALSE)
    p.cor<-as.numeric(unlist(p.cor)[2])
    bfc=alpha/length(p.cor)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    cex.cor.r <-cex.cor * abs(r)
    if (p.cor<=bfc){text(0.5, 0.5, txt, cex = 1.5,col='dark red')}else {text(0.5, 0.5, txt, cex = 1)}
}

#Holme Correction
panel.hccor <- function(x, y, digits=2, prefix="", cex.cor,method='spearman',alpha=0.05)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,method=method)
    p.cor<-cor.test(x,y,method=method,alternative='two.sided',exact=FALSE)
    p.cor<-as.numeric(unlist(p.cor)[2])
    p.cor=p.adjust(p.cor,method='holm')
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    txt <- paste(txt,round(p.cor,3),sep=', ')
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    cex.cor.r <-cex.cor * abs(r)
    if (p.cor<=alpha){text(0.5, 0.5, txt, cex = 1.5,col='dark red')}else {text(0.5, 0.5, txt, cex = 1)}
}


### put least square regression lines in the panel
#black least squares lines
panel.lm <- function (x, y, ...) {
	points(x, y, pch=20, ...) # this is what the default panel function does (except for pch=)
	abline(lm(y~x), lwd=2, col='black') # trend line from LOWESS smoother
}

#white least squares lines
panel.lmw <- function (x, y, ...) {
	points(x, y, pch=20, ...) # this is what the default panel function does (except for pch=)
	abline(lm(y~x), lwd=2, col="white") # trend line from LOWESS smoother
}

#red least squares lines
panel.lmr <- function (x, y, ...) {
	points(x, y, pch=20, ...) # this is what the default panel function does (except for pch=)
	abline(lm(y~x), lwd=2, col='red') # trend line from LOWESS smoother
}


#Plot the R-square values colored and scaled by a cutoff of 0.25
panel.r2 <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)^2
    r.col<-rep(2,length(r))
    r.col[r<0.25]<-1
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
#   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#    text(0.5, 0.5, txt, cex = cex.cor * r )
    text(0.5, 0.5, txt, cex = r.col*0.85,col=r.col)
}

#panel.cor2
#cor with cex=1.5
panel.cor2 <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
#    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#    text(0.5, 0.5, txt, cex = cex.cor * r )
    text(0.5, 0.5, txt, cex = 1.5)
}