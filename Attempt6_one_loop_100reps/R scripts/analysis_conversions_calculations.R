#analysis_conversions_calculations
#this script contains the conversions for files from Arlequin format to genepop format
#and loops that analyze the data involving converting genepop files to genind objects, 
#calculating the proportion of alleles captured, expected heterozygosity, and total number of alleles

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#defining root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt6_one_loop_100reps\\Simulations"
setwd(my_dir)

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig_highSamp", "\\lowMig_highSamp", "\\highMig_lowSamp", "\\lowMig_lowSamp")

#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1",
              "\\scen2",
              "\\scen3",
              "\\scen4",
              "\\scen5",
              "\\scen6",
              "\\scen7",
              "\\scen8",
              "\\scen9")

#import function
#converts arlequin files to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

#loop converting .arp to .gen for all combinations and scenarios
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    import_arp2gen_files(paste(my_dir,combinations[i],scenarios[j],sep=""), ".arp$")
  }
}

#DEFINING ARRAYS TO STORE RESULTS
#stores proportion of alleles captured
#high intensity
results_highMig_highSamp_equal = array(0, dim = c(9,100))
results_lowMig_highSamp_equal = array(0, dim = c(9,100))
results_highMig_highSamp_prop = array(0, dim = c(9,100))
results_lowMig_highSamp_prop = array(0, dim = c(9,100))
#low intensity
results_highMig_lowSamp_equal = array(0, dim = c(9,100))
results_lowMig_lowSamp_equal = array(0, dim = c(9,100))
results_highMig_lowSamp_prop = array(0, dim = c(9,100))
results_lowMig_lowSamp_prop = array(0, dim = c(9,100))

#stores the total alleles captured
total_alleles_highMig = array(0, dim = c(9, 100))
total_alleles_lowMig = array(0, dim = c(9, 100))
#defining array to store the expected heterozygosity for each locus
hexp_highMig = array(0, dim = c(9,100))
hexp_lowMig = array(0, dim = c(9,100))

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  #loop going through each scenario for each combination
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #loop going through every replicate for each scenario
    for(k in 1:length(list_files)) {
      #creating a temporary genind object
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining population boundaries by the first individual and the last individuals
      last_ind = as.numeric(cumsum(table(temp_genind@pop)))
      first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
      first_ind = first_ind[1:5]
      
      #defining how many individuals to sample from every population
      #high intensity
      if((i == 1) || (i == 2)) {
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.1)
        sample_size_equal = c(30,30,30,30,30)
      } else { #low intensity
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.5)
        sample_size_equal = c(15,15,15,15,15)
      }
      
      #defining individuals to randomly sample from each population for both strategies
      for(x in 1:5) {
        rows_to_samp_equal = sample(first_ind[x]:last_ind[x], sample_size_equal[x])
        rows_to_samp_prop = sample(first_ind[x]:last_ind[x], sample_size_prop[x])
      }
      
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
      
      #saving the total alleles present
      total_alleles = ncol(temp_genind@tab)
      
      #saving results 
      if(i == 1) {
        #saving proportion of alleles captured for both equal and proportional strategies
        results_highMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_highMig[j,k] = total_alleles
        
        #saving expected heterozygosity
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 2) { #if low migration
        #saving proportion of alleles captured 
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_lowMig[j,k] = total_alleles
        
        #saving expected heterozyogisty
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 3) {
        #saving prop alleles
        results_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
      } else {
        #saving prop alleles
        results_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
      }
    }
  }
}

#saving results
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt6_one_loop_100reps\\R scripts")
save(results_highMig_highSamp_equal, results_highMig_highSamp_prop, file="results_highMig_highSamp.Rdata")
save(results_lowMig_highSamp_equal, results_lowMig_highSamp_prop, file="results_lowMig_highSamp.Rdata")
save(results_highMig_lowSamp_equal, results_highMig_lowSamp_prop, file="results_highMig_lowSamp.Rdata")
save(results_lowMig_lowSamp_equal, results_lowMig_lowSamp_prop, file="results_lowMig_lowSamp.Rdata")

