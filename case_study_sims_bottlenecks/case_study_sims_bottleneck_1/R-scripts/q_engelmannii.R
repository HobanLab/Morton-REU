#Code for Q. engelmannii case study
#Code written by Kaylee Rosenberger, Emily Schumacher, and Dr. Sean Hoban in collaboration
#This script first imports simulation data and converts Arlequin file format to genepop format
#Next, the script loops over simulation replicates, converts each one into a temporary genind object (using Adegenet package),
#and simulates sampling from the wild populations by selecting a random number of individuals 
#The results are saved for each replicate in a matrix, run through a series of conversions/edits, and plotted using ggplot2

#Note: for these simulations, we were only interested in low sampling intensity
#however, increasing the intensity is as simple as increasing the values indicated on lines 58 and 59

#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)

#Flags 
#File conversion flag
#This flag is set to true when simulations have been run and files have been converted already
#There is no need to re-convert the files once they have been converted once
#if you want to re-run conversions, set this to FALSE
imported = FALSE

#set working directory
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims_bottlenecks\\case_study_sims_bottleneck_1\\Simulations\\q_engelmannii'
setwd(mydir)

#Defining an import function
#converts all arlequin files in a folder to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

#converting all simulation files from arlequin format to genepop format using defined import function
if(imported == FALSE) {
  import_arp2gen_files(mydir, ".arp$")
}

#pre-defining arrays to store results in loops below
#arrays dimensions are [1x100]
#the 100 represents simulation replicates
#the 1 represents 1 case study species with no varying parameters
#high intensity
results_q_engelmannii_equal_high = array(0, dim = c(1, 100))
results_q_engelmannii_prop_high = array(0, dim = c(1, 100))
#low intensity 
results_q_engelmannii_equal_low = array(0, dim = c(1, 100))
results_q_engelmannii_prop_low = array(0, dim = c(1, 100))

#high intensity
alleles_cap_q_engel_equal_high = array(0, dim = c(1, 100))
alleles_cap_q_engel_prop_high = array(0, dim = c(1, 100))
#low intensity
alleles_cap_q_engel_equal_low = array(0, dim = c(1, 100))
alleles_cap_q_engel_prop_low = array(0, dim = c(1, 100))

#array to store total alleles for each replicate
total_alleles_q_engelmannii = array(0, dim = c(1, 100))

#saving a list of the genind objects
quen_genind_list <- list() 

#list of hierfstat
quen_hierfstat <- list()

##quen pwfst array
quen_pwfst_array <- array(dim = c(4,4,100))

##min, max, mean of replicates 
quen_mean_max_min_fst <- matrix(nrow = 3, ncol = 100)

#***********************************************************************
#Flag
#flag set to true when you want to run fst
#set to false if you don't want it to run
f = FALSE

#Loop to simulate sampling
#First, create a list of all genepop files (all replicates) to loop over
#the variable 'i' represents each replicate
list_files = list.files(mydir, pattern = ".gen$")
for(x in 1:2) {
  for(i in 1:length(list_files)) {
    #creating a temporary genind object (using Adegenet package) for each simulation replicate
    temp_genind = read.genepop(list_files[[i]], ncode=3)
    
    if(f == TRUE) {
      ##creating genind list for QUEN genind 
      quen_genind_list[[i]] <- temp_genind
    
      ##convert genind files to hierfstat format to run pwfst 
      quen_hierfstat[[i]] <- genind2hierfstat(quen_genind_list[[i]])
    
      ##calculate statistics for QUEN - max, min, mean fst 
      quen_pwfst_array[,,i] <- pairwise.neifst(quen_hierfstat[[i]])
    
      ##calculate statistics for QUEN
      quen_mean_max_min_fst[1,i] <- mean(quen_pwfst_array[,,i], na.rm = TRUE)
      quen_mean_max_min_fst[2,i] <- min(quen_pwfst_array[,,i], na.rm = TRUE)
      quen_mean_max_min_fst[3,i] <- max(quen_pwfst_array[,,i], na.rm = TRUE)
    }
    
    #defining population boundaries by the first individual and the last individuals in each population
    #last individual for every population as the cumulative sum of all populations (ie., last individual for pop 1 is the sum of pop 1)
    last_ind = as.numeric(cumsum(table(temp_genind@pop)))
    #first individual of every population begins at 1, then for following populations, it is the last individual (cumulative sum) + 1
    #for example, if the last individual for pop 1 is 30, the first individual for pop 2 would be 31
    first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
    ##selecting the first 4 values since we have 4 populations
    first_ind = first_ind[1:4]
    
    #defining sample sizes for equal sampling and proportional sampling
    if(x==1) { #high intensity sampling
      sample_size_equal = c(100,100,100,100) #low intensity - if high, double values
      sample_size_prop = as.numeric(table(temp_genind@pop)*0.02)#low intensity - if high, double values
      sample_size_prop = ceiling(sample_size_prop) #round up decimal values to whole numbers
    }
    if(x==2) { #low intensity sampling
      sample_size_equal = c(50,50,50,50) #low intensity
      sample_size_prop = as.numeric(table(temp_genind@pop)*0.01)#low intensity
      sample_size_prop = ceiling(sample_size_prop) #round up decimal values to whole numbers
    }
    
    #defining 'rows' or individuals to sample from 
    rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]))
    rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]))
    
    #calculating and saving the raw alleles captured by each strategy into a variable
    sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
    sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
    
    #calculating the total alleles present
    total_alleles = ncol(temp_genind@tab)
    
    #saving the proportion of alleles captured by each strategy
    #ie., alleles sampled/total = proportion of alleles captured
    if(x==1) {
      results_q_engelmannii_equal_high[1,i] = sample_n_alleles_equal/total_alleles
      results_q_engelmannii_prop_high[1,i] = sample_n_alleles_prop/total_alleles
      
      alleles_cap_q_engel_equal_high[1,i] = sample_n_alleles_equal
      alleles_cap_q_engel_prop_high[1,i] = sample_n_alleles_prop
    }
    if(x==2) {
      results_q_engelmannii_equal_low[1,i] = sample_n_alleles_equal/total_alleles
      results_q_engelmannii_prop_low[1,i] = sample_n_alleles_prop/total_alleles
      
      alleles_cap_q_engel_equal_low[1,i] = sample_n_alleles_equal
      alleles_cap_q_engel_prop_low[1,i] = sample_n_alleles_prop
    }
     #saving the total alleles present
    total_alleles_q_engelmannii[1,i] = total_alleles
    
  }
}
#**************************************************************************************
#Preparing data for graphics
#converting results from matrices to data frames 
results_q_engelmannii_equal_high = as.data.frame(results_q_engelmannii_equal_high)
results_q_engelmannii_prop_high = as.data.frame(results_q_engelmannii_prop_high)

results_q_engelmannii_equal_low = as.data.frame(results_q_engelmannii_equal_low)
results_q_engelmannii_prop_low = as.data.frame(results_q_engelmannii_prop_low)

#converting data to long format
results_q_engelmannii_equal_high = gather(results_q_engelmannii_equal_high, replicate, prop_all)
results_q_engelmannii_prop_high = gather(results_q_engelmannii_prop_high, replicate, prop_all)

results_q_engelmannii_equal_low = gather(results_q_engelmannii_equal_low, replicate, prop_all)
results_q_engelmannii_prop_low = gather(results_q_engelmannii_prop_low, replicate, prop_all)

#variables to keep track of the strategy used -- equal or proportional
equal_strategy = rep("equal", 100)
prop_strategy = rep("proportional", 100)

#defining a column to keep track of the strategy variable
results_q_engelmannii_equal_high$strategy = equal_strategy
results_q_engelmannii_prop_high$strategy = prop_strategy

results_q_engelmannii_equal_low$strategy = equal_strategy
results_q_engelmannii_prop_low$strategy = prop_strategy

#combining both dataframes (holding equal and proportional data) into one large dataframe for plotting
combined_q_engelmannii_high = rbind(results_q_engelmannii_equal_high, results_q_engelmannii_prop_high)
combined_q_engelmannii_low = rbind(results_q_engelmannii_equal_low, results_q_engelmannii_prop_low)

#saving data to Rdata files
setwd('C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims_bottlenecks\\case_study_sims_bottleneck_1\\R-scripts')
save(results_q_engelmannii_equal_high, results_q_engelmannii_prop_high, results_q_engelmannii_equal_low, results_q_engelmannii_prop_low, file="results_q_engelmannii.Rdata")
save(combined_q_engelmannii_high, combined_q_engelmannii_low, file="combined_results_q_engelmannii.Rdata")

save(alleles_cap_q_engel_equal_high, alleles_cap_q_engel_prop_high, alleles_cap_q_engel_equal_low, alleles_cap_q_engel_prop_low, file="alleles_capt_q_engel.Rdata")

if(f == TRUE) {
save(quen_mean_max_min_fst, file="q_engelmannii_fst.Rdata")
}
#****************************************************************************************
#Creating Q. engelmannii graph using ggplot 2
#Note: this is only the results of the Q. engelmannii case study
#high intensity
p = ggplot(combined_q_engelmannii_high, aes(x=strategy, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ggtitle("Q. engelmannii") +
  xlab("Strategy") +
  ylab("Proportion of alleles captured") +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw() +
  theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))

#low intensity
p = ggplot(combined_q_engelmannii_low, aes(x=strategy, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ggtitle("Q. engelmannii") +
  xlab("Strategy") +
  ylab("Proportion of alleles captured") +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw() +
  theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))
