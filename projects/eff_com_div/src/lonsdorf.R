###Lonsdorf model
###Matlab code for simulations described in Shuster et al. Community Heritability MS

reps <- 10 #number of times to run simulation
### tree genotypic values - range from 11 to 21
trees2 <- read.csv('../data/lonsdorf/trees.txt') #SEE BELOW CODE FOR VALUES
###generate arthropod alleles for 25 species
insects <- matrix(NA,nrow=25,ncol=2)
insects[,1] <- runif(25,5,21) #heterozygote value between 5 and 21
insects[,2] <- runif(25,0,3) #range between 0 and 3
insect <- insects
insect[,1] <- (insects[,1]-0.5*insects[,2])/2 #C allelic value = (HET - 0.5*range)/2
insect[,2] <- (insects[,1]+0.5*insects[,2])/2 #D allelic value = (HET + 0.5*range)/2
###Experimental design
GG <- 8 #number of selection scenarios
YY <- 5 #number of environmental scenarios
###Parameter values
VeT <- 8 #environmental variance in tree trait influences tree heritability (2 for high H2 and 8 for low H2)
Ve <- 0.1 #environmental variance for insect trait
K <- 100 #insect pop carrying capacity
VeN <- 15 #step size for environmental variance in interactions (0 15 30 45 60)
T <- nrow(trees2) #number of trees
I <- nrow(insect) #number of insects
art_g <- matrix(NA,nrow=T,ncol=I) #insect growth rate
art_z <- matrix(NA,nrow=T,ncol=I) #insect ?
art_Vg <- matrix(NA,nrow=T,ncol=I) #insect ?
art_Vz <- matrix(NA,nrow=T,ncol=I) #insect ?
gen_load <- matrix(NA,nrow=T,ncol=I) #insect ?
dem_load <- matrix(NA,nrow=T,ncol=I) #insect ?
art_pop <- matrix(NA,nrow=T,ncol=I) #insect ?
                                        #STORAGE LIST
out <- gg.list <- yy.list <- list()
###Simulation
tic <- Sys.time() #start stopwatch
for (RR in 1:reps){
###generating trees to use
  scores_XX <- matrix(NA,nrow=nrow(trees2),ncol=GG)
  scores_XX[,1] <- trees2[,2] + runif(T,0,1) * VeT - VeT/2 #scores are in second column of trees
###filling phenotypic values for trees (same set of trees for all scenarios
  for (z in 2:GG){
    scores_XX[1:T,z] <- scores_XX[1:T,1]
  }
  trees <- scores_XX #just renaming phenotypes to be trees
  for (y in 1:YY){
                                        #YY variation scenarios of other ecological interactions
    for (z in 1:GG){
                                        #GG selection intensity scenarios
      for (i in 1:T){
                                        #insects on trees
                                        #for each tree i
        for (j in 1:I){
                                        #for each insect j
###Equation 6 from MS - Arthropod (art) Gene frequency
          if (trees[i,z] < 2*insect[j,1]){
            art_g[i,j] <- 0
          }else if (trees[i,z] > 2*insect[j,2]){
            art_g[i,j] <- 1
          }else{
            art_g[i,j] <- (trees[i,z] - 2*insect[j,1]) / (2*insect[j,2] - 2*insect[j,1])
          }
###Equation 6 from MS - calculating mean trait Z from art gene frequency
          art_z[i,j] <- 2*insect[j,2]*art_g[i,j]^2 + 2*art_g[i,j]*(1-art_g[i,j])*(insect[j,1]+insect[j,2])+2*insect[j,1]*(1-art_g[i,j])^2 + runif(1)*Ve - Ve/2
###Art genetic (Vg) and trait variance (Vz)
          art_Vg[i,j] <- 2*art_g[i,j]*(1-art_g[i,j])
          art_Vz[i,j] <- art_Vg[i,j]*(insect[j,2] - insect[j,1])^2+Ve
###Evolutionary (gen) and demographic (dem) loads from selection
          gen_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_z[i,j]- trees[i,z])^2
          dem_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_Vz[i,j])
###Equation 7 from MS - art predicted population size as a function of loads and ecological variance
          art_pop[i,j] <- K * (1 - gen_load[i,j] - dem_load[i,j])+runif(1)*VeN*(y-1)-VeN*(y-1)/2
###preventing art pops from going negative or zero slightly above zero
          if (art_pop[i,j] < 0){
            art_pop[i,j] <- runif(1)*3;
          }else{
            art_pop[i,j] #not sure if this is correct
          }
        } #end insect loop
      } #end tree loop
###keeping track of which iteration the simulation is on reps*GG*(y-1)+reps*(z-1)+RR
###dlmwrite(int2str(reps*GG*(y-1)+reps*(z-1)+RR),art_pop,'\t');  %%% write arthropod sampling to file %%%
      print(paste(reps*GG*(y-1),reps*(z-1),RR,sep='_'))
      gg.list[[z]] <- art_pop
    } #end selection loop
    yy.list[[y]] <- gg.list
  } #end environment loop
  toc <- Sys.time() #stop stopwatch
  out[[RR]] <- yy.list
} #end REP loop and END simulation
###4.344779 per rep
print(toc-tic)
dput(out,file='./ld.Rdata')
