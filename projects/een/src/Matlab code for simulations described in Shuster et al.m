%% Matlab code for simulations described in Shuster et al.’s Community Heritability MS

%% Comments are designated by “%”

clear all
close all

%%% tree genotypic values – range from 11 to 21 %%%
trees2=csvread('../data/trees.txt'); %%% SEE BELOW CODE FOR VALUES
trees2=trees2(2:51,:);


%%% generate arthropod alleles for 25 species %%%
insects(:,1)=5+rand(25,1)*16;  %%% heterozygote value between 5 and 21 
insects(:,2)=rand(25,1)*3;     %%% range between 0 and 3
insect(:,1)=(insects(:,1)-0.5*insects(:,2))/2; %%% C allelic value = (HET - 0.5*range)/2
insect(:,2)=(insects(:,1)+0.5*insects(:,2))/2; %%% D allelic value = (HET + 0.5*range)/2

%%% Experimental design %%%
reps = 10; %% number of times to run simulation
GG = 8; %% number of selection scenarios
YY = 5; %% number of environmental scenarios

%%% parameter values %%%
VeT = 8;  %% environmental variance in tree trait influences tree heritability (2 for high H2 and 8 for low H2)
Ve = .1; %% environmental variance for insect trait
K = 100;  %% insect pop carrying capacity
VeN = 15; %% step size for environmental variance in interactions (0 15 30 45 60)

art_r = 1; %% insect growth rate
T=size(trees2,1);
I=size(insect,1);


%%% simulation begins %%%
for RR = 1:reps

  tic %% start stop watch

%%% generating trees to use %%%
  scores_XX(:,1) = trees2(:,2)+rand(T,1)*VeT-VeT/2; %%% scores are in second column of trees

%%% filling phenotypic values for trees (same set of trees for all
%%% scenarios
  for z=2:GG;
    scores_XX(1:T,z) = scores_XX(1:T,1);    
  end
  trees = scores_XX(:,:); %% just renaming phenotypes to be trees


  for y =1:YY  %% YY VARIATION SCENARIOS OF OTHER ECOLOGICAL INTERACTIONS

    for z=1:GG;  %% GG SELECTION INTENSITY SCENARIOS

      %% INSECTS ON TREES
      for i=1:T  %for each tree i
        
	for j = 1:I  %% for each insect j
          
%%% Equation 6 from MS - Arthropod (art) Gene frequency  %%%
	  if trees(i,z) < 2*insect(j,1)
	    art_g(i,j) = 0;
	  elseif trees(i,z) > 2*insect(j,2)
	    art_g(i,j) = 1;
	  else art_g(i,j) = (trees(i,z) - 2*insect(j,1))/(2*insect(j,2) - 2*insect(j,1));
          end
          
%%% Equation 6 from MS - calculating mean trait Z from art gene frequency %%%
	  art_z(i,j) = 2*insect(j,2)*art_g(i,j)^2 + 2*art_g(i,j)*(1-art_g(i,j))*(insect(j,1)+insect(j,2))+2*insect(j,1)*(1-art_g(i,j))^2 + rand*Ve - Ve/2; 
          
%%% Art genetic (Vg) and trait variance (Vz) %%%
	  art_Vg(i,j) = 2*art_g(i,j)*(1-art_g(i,j));
	  art_Vz(i,j) = art_Vg(i,j)*(insect(j,2) - insect(j,1))^2+Ve;
          
%%% evolutionary (gen) and demographic (dem) loads from selection %%%
	  gen_load(i,j)=0.5*(0.00007924*2.511886^(z-1))*(art_z(i,j)- trees(i,z))^2;
	  dem_load(i,j)=0.5*(0.00007924*2.511886^(z-1))*(art_Vz(i,j));
          
%%% Equation 7 from MS - art predicted population size as a function of loads
%%% and ecological variance
	  art_pop(i,j) = K * (1 - gen_load(i,j) - dem_load(i,j))+rand*VeN*(y-1)-VeN*(y-1)/2;
          
%%% preventing art pops from going negative or zero
%%% (slightly above zero)
	  if art_pop(i,j) < 0
	    art_pop(i,j) = rand*3; 
	  else art_pop(i,j) ;
          end
          
        end  %% end insect loop
      end  %% end tree loop
      
%%% keeping track of which iteration the simulation is on %%%
      reps*GG*(y-1)+reps*(z-1)+RR
      k = k+1
      dlmwrite(strcat('../data/lonsdorf_out/',int2str(k),art_pop,'\t');  %%% write arthropod sampling to file %%%   
%%      dlmwrite(strcat('../data/lonsdorf_out/',int2str(reps*GG*(y-1)+reps*(z-1)+RR)),art_pop,'\t');  %%% write arthropod sampling to file %%%   
      
    end  %% end selection loop

  end  %% end environment loop

  toc %% stop stopwatch

end %%% end rep loop - END SIMULATION


Tree.txt table used in the code above
First column – clone ID; Second column – genotypic value (50 trees, 10 genotypes x 5 reps per genotype)

1 11
1 11
1 11
1 11
1 11
2 12.5
2 12.5
2 12.5
2 12.5
2 12.5
3 13.75
3 13.75
3 13.75
3 13.75
3 13.75
4 16
4 16
4 16
4 16
4 16
5 14
5 14
5 14
5 14
5 14
6 15.25
6 15.25
6 15.25
6 15.25
6 15.25
7 17.5
7 17.5
7 17.5
7 17.5
7 17.5
8 16.5
8 16.5
8 16.5
8 16.5
8 16.5
9 18.75
9 18.75
9 18.75
9 18.75
9 18.75
10 21
10 21
10 21
10 21
10 21
