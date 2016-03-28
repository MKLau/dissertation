source('../src/loadPitdata.R')
out.loc <- tail(commandArgs(),2)
x <- tree.nets[[as.numeric(args[2])]]
diag(x) <- 0
out <- computeModules(x)
out <- slot(out,"likelihood")
write.table(out.loc,append=T)
