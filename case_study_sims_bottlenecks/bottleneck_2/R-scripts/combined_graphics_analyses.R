#Combined_graphics.R
#Code written by Kaylee Rosenberger, Emily Schumacher, and Dr. Sean Hoban in collaboration
#This script combines all results from three case study species
#into one large dataframe, so that results could be plotted on one graph for comparison.

#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)

#loading in data that was saved from other case study R scripts
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims_bottlenecks\\bottleneck_2\\R-scripts'
setwd(mydir)
load("combined_results_q_acerifolia.Rdata")
load("q_acerifolia_fst.Rdata")

load("combined_results_q_engelmannii.Rdata")
load("q_engelmannii_fst.Rdata")

load("combined_results_q_oglethorpensis.Rdata")
load("q_oglethorpensis_fst.Rdata")

###LOW INTENSITY SAMPLING PLOT
#preparing data by naming a column indicating each case study species
#repeating 200 times for each species since there are 100 simulation replicates that were analzyed for equal and proportional (200 entries total)
q_acer = rep("Q. acerifolia", 200)
combined_q_acerifolia_low$species = q_acer

q_engel = rep("Q. engelmannii", 200)
combined_q_engelmannii_low$species = q_engel

q_ogle = rep("Q. oglethorpensis", 200)
combined_q_oglethorpensis_low$species = q_ogle

#combining all results into one larger dataframe for plotting/comparison on one graph
all_case_studies_low = rbind(combined_q_acerifolia_low, combined_q_engelmannii_low, combined_q_oglethorpensis_low)

#Creating graph using ggplot2 with all case study species compared together
#x-axis: species, y-axis: proportion of alleles captured
#the color of the boxplot indicates which strategy was used
p = ggplot(all_case_studies_low, aes(x=species, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) + #this line shows stars *** for significant values on the plot
  ggtitle("Case study species (Bottleneck 2)") + #labels for plot
  xlab("Species") +
  ylab("Proportion of alleles captured") +
  ylim(0.7,1.02) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) + #design elements
  theme_bw() #+
  #theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14)) #creating/displaying the plot and changing font size to be larger

###HIGH INTENSITY SAMPLING PLOT
#preparing data by naming a column indicating each case study species
#repeating 200 times for each species since there are 100 simulation replicates that were analzyed for equal and proportional (200 entries total)
q_acer = rep("Q. acerifolia", 200)
combined_q_acerifolia_high$species = q_acer

q_engel = rep("Q. engelmannii", 200)
combined_q_engelmannii_high$species = q_engel

q_ogle = rep("Q. oglethorpensis", 200)
combined_q_oglethorpensis_high$species = q_ogle

#combining all results into one larger dataframe for plotting/comparison on one graph
all_case_studies_high = rbind(combined_q_acerifolia_high, combined_q_engelmannii_high, combined_q_oglethorpensis_high)

#Creating graph using ggplot2 with all case study species compared together
#x-axis: species, y-axis: proportion of alleles captured
#the color of the boxplot indicates which strategy was used
p = ggplot(all_case_studies_high, aes(x=species, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) + #this line shows stars *** for significant values on the plot
  ggtitle("Case study species (Bottleneck 2)") + #labels for plot
  xlab("Species") +
  ylab("Proportion of alleles captured") +
  ylim(0.7,1.02) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) + #design elements
  theme_bw() #+
#theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14)) #creating/displaying the plot and changing font size to be larger


###################################################################################
#Calculating Fst of each species 
#only run if f == TRUE for all case study species
##final df with combined species
species_pwfst_df <- matrix(nrow = 3, ncol = 3)
###write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
for(a in 1:length(quog_mean_max_min_fst[,1])){
  species_pwfst_df[1,a] <-  round(mean(quac_mean_max_min_fst[a,]),3)  
  species_pwfst_df[2,a] <-  round(mean(quen_mean_max_min_fst[a,]),3)
  species_pwfst_df[3,a] <- round(mean(quog_mean_max_min_fst[a,]),3)
}

##rename dfs 
rownames(species_pwfst_df) <- c("QUAC","QUEN","QUOG")
colnames(species_pwfst_df) <- c("MeanFst", "MinFst", "MaxFst")

##write out to csv
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims_bottlenecks\\bottleneck_2\\R scripts'
setwd(mydir)
write.csv(species_pwfst_df, "case_study_fst.csv")
