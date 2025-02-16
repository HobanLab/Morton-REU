#analysis_conversions_calculations R script
#Code written by Kaylee Rosenberger, Emily Schumacher, and Dr. Sean Hoban in collaboration
#this script contains the conversions for files from Arlequin format to genepop format
#and loops that analyze the data, involving converting genepop files to genind objects, 
#calculating the proportion of alleles captured, expected heterozygosity, and total number of alleles
#and saving those results into matrices

#Library functions
library(adegenet)
library(hierfstat)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)

#Flags
#File conversion flag
#This flag is set to false when simulations have been run and files have been converted already
#There is no need to re-convert the files once they have been converted once
#if you want to re-run conversions or you re-ran simulations, set this to TRUE
convert = FALSE 
#Fst flag
#Fst code adds a lot of time to run the code 
#so if you don't want to run it, keep Fst off by setting it FALSE
f <- FALSE

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
if(convert == TRUE) {
  for(i in 1:length(combinations)) {
    for(j in 1:length(scenarios)) {
      import_arp2gen_files(paste(my_dir,combinations[i],scenarios[j],sep=""), ".arp$")
    }
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

#stores the number of alleles captured by the sample strategy (instead of proportion of alleles captured)
#high intensity
samp_all_highMig_highSamp_equal = array(0, dim = c(9,100))
samp_all_lowMig_highSamp_equal = array(0, dim = c(9,100))
samp_all_highMig_highSamp_prop = array(0, dim = c(9,100))
samp_all_lowMig_highSamp_prop = array(0, dim = c(9,100))
#low intensity
samp_all_highMig_lowSamp_equal = array(0, dim = c(9,100))
samp_all_lowMig_lowSamp_equal = array(0, dim = c(9,100))
samp_all_highMig_lowSamp_prop = array(0, dim = c(9,100))
samp_all_lowMig_lowSamp_prop = array(0, dim = c(9,100))

#defining array to store the expected heterozygosity for each locus
hexp_highMig = array(0, dim = c(9,100))
hexp_lowMig = array(0, dim = c(9,100))

#####Fst saving
##high mig lists
highmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                          list(), list(), list(), list())

##low mig lists
lowmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                         list(), list(), list(), list())

##arrays low and high mig 
highmig_pwfst_array <- array(dim = c(5,5,100,9))
lowmig_pwfst_array <- array(dim = c(5,5,100,9))

##min max mean fst arrays
highmig_fst_min_mean_max <- array(dim = c(3,100,9))
lowmig_fst_min_mean_max <- array(dim = c(3,100,9))

##results dfs
highmig_pwfst_output <- matrix(nrow = 9, ncol = 3)
lowmig_pwfst_output <- matrix(nrow = 9, ncol = 3)

###########################################################################################################
#Sampling loop - converts replicates to genind objects and performs sampling 
#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
#first loop over combinations (migration rate + sample intensity)
#i = combination
for(i in 1:length(combinations)) {
  #next, loop through each scenario (1-9) for each combination
  #j = scenario
  for(j in 1:length(scenarios)) {
    #update working directory for each scenario (to navigate through files)
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    #create a list of all .gen files (ie., simulation replicates)
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #finally, loop through every replicate for each scenario, in each combination
    #k = replicate
    for(k in 1:length(list_files)) {
      #creating a temporary genind object of the simulation replicate
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining population boundaries by the first individual and the last individuals
      #last individual for every population as the cumulative sum of all populations (ie., last individual for pop 1 is the sum of pop 1)
      last_ind = as.numeric(cumsum(table(temp_genind@pop)))
      #first individual of every population begins at 1, then for following populations, it is the last individual (cumulative sum) + 1
      #for example, if the last individual for pop 1 is 30, the first individual for pop 2 would be 31
      first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
      #selecting the first 5 values since we have 5 populations
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
      #from above, we defined the population boundaries of each population by the first and last individuals of the populations
      #ie., the 'range' of the populations is defined as (first_ind:last_ind)
      #here, we are defining which individuals to sample from, for each of the populations
      #so we are sampling a certain amount (define in sample_size... variable) from the population 'ranges'
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
        
        #saving number of alleles captured
        samp_all_highMig_highSamp_equal[j,k] = sample_n_alleles_equal
        samp_all_highMig_highSamp_prop[j,k] = sample_n_alleles_prop
        
        #saving total alleles present for each replicate
        total_alleles_highMig[j,k] = total_alleles
        
        #saving expected heterozygosity
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
       
        if(f == TRUE) {
        ##convert to hierfstat data frame and store in list
        highmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array - 
        highmig_pwfst_array[,,k,j] <-  pairwise.neifst(highmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate
        highmig_fst_min_mean_max[1,k,j] <- mean(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[2,k,j] <- min(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[3,k,j] <- max(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
        for(a in 1:length(scenarios)){
          for(b in 1:3){
            highmig_pwfst_output[a,b] <- round(mean(highmig_fst_min_mean_max[b,,a]),3)
          }
        }
      }
      } else if (i == 2) { #if the loop is in the second combination (low migration, high intensity), save results here
        #saving proportion of alleles captured for both equal and proportional strategies
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving number of alleles captured
        samp_all_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal
        samp_all_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop
        
        #saving total alleles present
        total_alleles_lowMig[j,k] = total_alleles
        
        #saving expected heterozyogisty
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
        
        if(f == TRUE) {
        ##convert to hierfstat data frame and store in a list
        lowmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array 
        lowmig_pwfst_array[,,k,j] <-  pairwise.neifst(lowmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate
        lowmig_fst_min_mean_max[1,k,j] <- mean(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[2,k,j] <- min(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[3,k,j] <- max(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
        for(a in 1:length(scenarios)){
          for(b in 1:3){
            lowmig_pwfst_output[a,b] <- round(mean(lowmig_fst_min_mean_max[b,,a]),3)
          }
        }
       }
      } else if (i == 3) { #if loop is in high migration, low intensity combination, save results here
        ##saving proportion of alleles captured for both equal and proportional strategies
        results_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving number of alleles captured
        samp_all_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal
        samp_all_highMig_lowSamp_prop[j,k] =  sample_n_alleles_prop
        
        #we don't need to save the total alleles present or hexp again, since we are using the same data as the high migration combination above (i = 1)
        
      } else { # (i==4) if loop is in low migration, low intensity combination, save results here
        #saving proportion of alleles captured for both equal and proportional strategies
        results_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving number of alleles captured 
        samp_all_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal
        samp_all_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop
        
        #we don't need to save the total alleles present or hexp again, since we are using the same data as the low migration combination above (i = 2)
      }
    }
  }
}

#######################################################################################################################
#saving results to R data file
#save proportion of alleles captured 
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims\\Simulations")
save(results_highMig_highSamp_equal, results_highMig_highSamp_prop, file="results_highMig_highSamp.Rdata")
save(results_lowMig_highSamp_equal, results_lowMig_highSamp_prop, file="results_lowMig_highSamp.Rdata")
save(results_highMig_lowSamp_equal, results_highMig_lowSamp_prop, file="results_highMig_lowSamp.Rdata")
save(results_lowMig_lowSamp_equal, results_lowMig_lowSamp_prop, file="results_lowMig_lowSamp.Rdata")

#save number of alleles captured
save(samp_all_highMig_highSamp_equal, samp_all_highMig_highSamp_prop, file="alleles_capt_highMig_highSamp.Rdata")
save(samp_all_highMig_lowSamp_equal, samp_all_highMig_lowSamp_prop, file="alleles_capt_highMig_lowSamp.Rdata")
save(samp_all_lowMig_highSamp_equal, samp_all_lowMig_highSamp_prop, file="alleles_capt_lowMig_highSamp.Rdata")
save(samp_all_lowMig_lowSamp_equal, samp_all_lowMig_lowSamp_prop, file="alleles_capt_lowMig_lowSamp.Rdata")

#save total alleles present
save(total_alleles_highMig, total_alleles_lowMig, file="total_alleles.Rdata")

# #save Fst files
# save(highmig_fst_min_mean_max, file="highMig_fst.Rdata")
# save(highmig_pwfst_output, file="highMig_fst_output.Rdata")
# save(lowmig_fst_min_mean_max, file="lowMig_fst.Rdata")
# save(lowmig_pwfst_output, file="lowMig_fst_output.Rdata")
