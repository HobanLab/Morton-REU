#Code for the Q. oglethorpensis case study
#This script first imports simulation data and converts Arlequin file format to genepop format
#Next, the script loops over simulation replicates, converts each one into a temporary genind object (using Adegenet package),
#and simulates sampling from the wild populations by selecting a random number of individuals 
#The results are saved for each replicate in a matrix, run through a series of conversions/edits, and plotted using ggplot2

#library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#set working directory
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_studies\\Simulations\\q_oglethorpensis'
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
import_arp2gen_files(mydir, ".arp$")

##pre-defining arrays to store results from loop below
results_q_oglethorpensis_equal = array(0, dim = c(1, 100))
results_q_oglethorpensis_prop = array(0, dim = c(1, 100))

#array to store total alleles for each replicate
total_alleles_q_oglethorpensis = array(0, dim = c(1, 100))

#***********************************************************************
#Loop to simulate sampling
#First, create a list of all genepop files (all replicates) to loop over
#the variable 'i' represents each replicate
list_files = list.files(mydir, pattern = ".gen$")
for(i in 1:length(list_files)) {
  #creating a temporary genind object (using Adegenet package) for each simulation replicate 
  temp_genind = read.genepop(list_files[[i]], ncode=3)
  
  #defining population boundaries by the first individual and the last individuals in each population
  last_ind = as.numeric(cumsum(table(temp_genind@pop)))
  first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
  first_ind = first_ind[1:5]
  
  #defining sample sizes for equal sampling and proportional sampling
  sample_size_equal = c(10,10,10,10,10) #low intensity - if high, double values
  sample_size_prop = as.numeric(table(temp_genind@pop)*0.05)#low intensity - if high, double values
  sample_size_prop = ceiling(sample_size_prop) #round up decimal values to whole numbers
  
  #defining 'rows' or individuals to sample from 
  rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]), sample(first_ind[5]:last_ind[5], sample_size_equal[5]))
  rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]), sample(first_ind[5]:last_ind[5], sample_size_prop[5]))
  
  #calculating and saving the raw alleles captured by each strategy into a variable
  sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
  sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
  
  #saving the total alleles present
  total_alleles = ncol(temp_genind@tab)
  
  #saving the proportion of alleles captured by each strategy
  results_q_oglethorpensis_equal[1,i] = sample_n_alleles_equal/total_alleles
  results_q_oglethorpensis_prop[1,i] = sample_n_alleles_prop/total_alleles
  
  #saving the total alleles present
  total_alleles_q_oglethorpensis[1,i] = total_alleles
  
}

#**************************************************************************************
#Preparing data for graphics
#converting results from matrices to data frames 
results_q_oglethorpensis_equal = as.data.frame(results_q_oglethorpensis_equal)
results_q_oglethorpensis_prop = as.data.frame(results_q_oglethorpensis_prop)

#converting to long format
results_q_oglethorpensis_equal = gather(results_q_oglethorpensis_equal, replicate, prop_all)
results_q_oglethorpensis_prop = gather(results_q_oglethorpensis_prop, replicate, prop_all)

#variables to keep track of the strategy used -- equal or proportional
equal_strategy = rep("equal", 100)
prop_strategy = rep("proportional", 100)

#defining a column to keep track of the strategy variable
results_q_oglethorpensis_equal$strategy = equal_strategy
results_q_oglethorpensis_prop$strategy = prop_strategy

#combining both dataframes (holding equal and proportional data) into one large dataframe for plotting
combined_q_oglethorpensis = rbind(results_q_oglethorpensis_equal, results_q_oglethorpensis_prop)

#saving data to Rdata files
setwd(mydir)
save(results_q_oglethorpensis_equal, results_q_oglethorpensis_prop, file="results_q_oglethorpensis.Rdata")
save(combined_q_oglethorpensis, file="combined_results_q_oglethorpensis.Rdata")

#****************************************************************************************
#Creating graph for Q. oglethorpensis
#Note: this graph only shows results for Q. oglethorpensis
p = ggplot(combined_q_oglethorpensis, aes(x=strategy, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ggtitle("Q. oglethorpensis") +
  xlab("Strategy") +
  ylab("Proportion of alleles captured") +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw() +
  theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))
