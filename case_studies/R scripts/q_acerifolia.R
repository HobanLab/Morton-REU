#Q. Acerifolia

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_studies\\Simulations\\q_acerifolia'
setwd(mydir)

#import function
#converts arlequin files to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

#converting all simulation files from arlequin format to genepop format
import_arp2gen_files(mydir, ".arp$")

#defining arrays to store results
results_q_acerifolia_equal = array(0, dim = c(1, 100))
results_q_acerifolia_prop = array(0, dim = c(1, 100))

#array to store total arrays for each replicate
total_alleles_q_acerifolia = array(0, dim = c(1, 100))

#***********************************************************************
#Loop to simulate sampling
list_files = list.files(mydir, pattern = ".gen$")
for(i in 1:length(list_files)) {
  temp_genind = read.genepop(list_files[[i]], ncode=3)
  
  #defining population boundaries by the first individual and the last individuals
  last_ind = as.numeric(cumsum(table(temp_genind@pop)))
  first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
  first_ind = first_ind[1:4]
  
  sample_size_equal = c(8,8,8,8)
  sample_size_prop = as.numeric(table(temp_genind@pop)*0.06)
  sample_size_prop = ceiling(sample_size_prop) #round up values
  
  rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]))
  rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]))
  
  sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
  sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
  
  #saving the total alleles present
  total_alleles = ncol(temp_genind@tab)
  
  results_q_acerifolia_equal[1,i] = sample_n_alleles_equal/total_alleles
  results_q_acerifolia_prop[1,i] = sample_n_alleles_prop/total_alleles
  
  total_alleles_q_acerifolia[1,i] = total_alleles
  
}

#**************************************************************************************
#Preparing data for graphics
#converting results to data frames 
results_q_acerifolia_equal = as.data.frame(results_q_acerifolia_equal)
results_q_acerifolia_prop = as.data.frame(results_q_acerifolia_prop)
#converting to long format
results_q_acerifolia_equal = gather(results_q_acerifolia_equal, replicate, prop_all)
results_q_acerifolia_prop = gather(results_q_acerifolia_prop, replicate, prop_all)

#variables to keep track of strategy
equal_strategy = rep("equal", 100)
prop_strategy = rep("proportional", 100)

results_q_acerifolia_equal$strategy = equal_strategy
results_q_acerifolia_prop$strategy = prop_strategy

#combining dataframes into one
combined_q_acerifolia = rbind(results_q_acerifolia_equal, results_q_acerifolia_prop)

setwd(mydir)
save(results_q_acerifolia_equal, results_q_acerifolia_prop, file="results_q_acerifolia.Rdata")
save(combined_q_acerifolia, file="combined_results_q_acerifolia.Rdata")

#****************************************************************************************
#Creating graphs
p = ggplot(combined_q_acerifolia, aes(x=strategy, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ggtitle("Q. acerifolia") +
  xlab("Strategy") +
  ylab("Proportion of alleles captured") +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw()
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))
