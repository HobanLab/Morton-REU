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

#loading in data that was saved from other case study R scripts
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_acerifolia'
setwd(mydir)
load("combined_results_q_acerifolia.Rdata")

mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_engelmannii'
setwd(mydir)
load("combined_results_q_engelmannii.Rdata")

mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_oglethorpensis'
setwd(mydir)
load("combined_results_q_oglethorpensis.Rdata")

#preparing data by naming a column indicating each case study species
#repeating 200 times for each species since there are 100 simulation replicates that were analzyed for equal and proportional (200 entries total)
q_acer = rep("Q. acerifolia", 200)
combined_q_acerifolia$species = q_acer

q_engel = rep("Q. engelmannii", 200)
combined_q_engelmannii$species = q_engel

q_ogle = rep("Q. oglethorpensis", 200)
combined_q_oglethorpensis$species = q_ogle

#combining all results into one larger dataframe for plotting/comparison on one graph
all_case_studies = rbind(combined_q_acerifolia, combined_q_engelmannii, combined_q_oglethorpensis)

#Creating graph using ggplot2 with all case study species compared together
#x-axis: species, y-axis: proportion of alleles captured
#the color of the boxplot indicates which strategy was used
p = ggplot(all_case_studies, aes(x=species, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) + #this line shows stars *** for significant values on the plot
  ggtitle("Case study species") + #labels for plot
  xlab("Species") +
  ylab("Proportion of alleles captured") +
  ylim(0.7,1.02) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) + #design elements
  theme_bw() #+
  #theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14)) #creating/displaying the plot and changing font size to be larger
