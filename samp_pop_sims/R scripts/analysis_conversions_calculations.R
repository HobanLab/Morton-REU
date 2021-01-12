#analysis_conversions_calculations R script
#this script contains the conversions for files from Arlequin format to genepop format
#and loops that analyze the data, involving converting genepop files to genind objects, 
#calculating the proportion of alleles captured, expected heterozygosity, and total number of alleles
#and saving those results into matrices

#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#defining root directory (containing sub-folders)
#and setting working directory
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims\\Simulations"
setwd(my_dir)

#list of combinations
#combination sub-folder directories extensions
combinations = c("\\highMig_highSamp", "\\lowMig_highSamp", "\\highMig_lowSamp", "\\lowMig_lowSamp")

#list of scenarios
#simulation sub-folder directories extensions
scenarios = c("\\scen1",
              "\\scen2",
              "\\scen3",
              "\\scen4",
              "\\scen5",
              "\\scen6",
              "\\scen7",
              "\\scen8",
              "\\scen9")

#defining an import function
#converts all arlequin files in a folder to genepop files
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
#stores proportion of alleles captured for both strategies -- equal and proportional 
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

#stores the total alleles present
total_alleles_highMig = array(0, dim = c(9, 100))
total_alleles_lowMig = array(0, dim = c(9, 100))
#defining array to store the expected heterozygosity for each locus
hexp_highMig = array(0, dim = c(9,100))
hexp_lowMig = array(0, dim = c(9,100))

###########################################################################################################
#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
#first loop over combinations (migration rate + sample intensity)
for(i in 1:length(combinations)) {
  #next, loop through each scenario (1-9) for each combination
  for(j in 1:length(scenarios)) {
    #update working directory for each scenario (to navigate through files)
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    #create a list of all .gen files (ie., simulation replicates)
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #finally, loop through every replicate for each scenario, in each combination
    for(k in 1:length(list_files)) {
      #creating a temporary genind object of the simulation replicate
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining population boundaries by the first individual and the last individuals
      last_ind = as.numeric(cumsum(table(temp_genind@pop)))
      first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
      first_ind = first_ind[1:5]
      
      #defining how many individuals to sample from every population
      #it varies based on which combination the loop is in (based on sample intensity)
      #the first two combinations are high intensity sampling, the last two are low intensity
      #if high intensity, sample 10% of each population, or 30 individuals equally
      if((i == 1) || (i == 2)) {
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.1)
        sample_size_equal = c(30,30,30,30,30)
      } else { #if low intensity, sample 5% of each population , or 15 individuals equally
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.05)
        sample_size_equal = c(15,15,15,15,15)
      }
      
      #defining individuals to randomly sample from each population for both strategies
      rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]), sample(first_ind[5]:last_ind[5], sample_size_equal[5]))
      rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]), sample(first_ind[5]:last_ind[5], sample_size_prop[5]))
      
      #calculating and saving the alleles sampled
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
      
      #saving the total alleles present
      total_alleles = ncol(temp_genind@tab)
      
      #saving results in matrices
      #where the results need to be saved depends on what iteration the loop is in
      #if the loop is the first combination (high intensity, high migration), save results here (checking based on variable 'i', the loop counter for the outermost loop)
      if(i == 1) {
        #saving proportion of alleles captured for both equal and proportional strategies
        results_highMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles present for each replicate
        total_alleles_highMig[j,k] = total_alleles
        
        #saving expected heterozygosity
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 2) { #if the loop is in the second combination (low migration, high intensity), save results here
        #saving proportion of alleles captured for both equal and proportional strategies
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles present
        total_alleles_lowMig[j,k] = total_alleles
        
        #saving expected heterozyogisty
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 3) { #if loop is in high migration, low intensity combination, save results here
        ##saving proportion of alleles captured for both equal and proportional strategies
        results_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #we don't need to save the total alleles present or hexp again, since we are using the same data as the high migration combination above (i = 1)
        
      } else { #if loop is in low migration, low intensity combination, save results here
        #saving proportion of alleles captured for both equal and proportional strategies
        results_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #we don't need to save the total alleles present or hexp again, since we are using the same data as the low migration combination above (i = 2)
      }
    }
  }
}

#######################################################################################################################
#saving results to R data file
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt6_one_loop_100reps\\R scripts")
save(results_highMig_highSamp_equal, results_highMig_highSamp_prop, file="results_highMig_highSamp.Rdata")
save(results_lowMig_highSamp_equal, results_lowMig_highSamp_prop, file="results_lowMig_highSamp.Rdata")
save(results_highMig_lowSamp_equal, results_highMig_lowSamp_prop, file="results_highMig_lowSamp.Rdata")
save(results_lowMig_lowSamp_equal, results_lowMig_lowSamp_prop, file="results_lowMig_lowSamp.Rdata")
