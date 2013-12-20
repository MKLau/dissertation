###Plot the trajectory of each arthropod on each tree over the course of the simulation
ldr <- dget('./ld.Rdata')
f <- read.csv('trees.txt')
                                        #REP - Environment - Selection
AP <- ldr[[1]][[1]]

library(abind)
AP. <- do.call(abind,c(AP,along=3))
                                        #mean abundances by tree
ma.tree <- apply(AP.,c(1,3),mean)
