#libraries
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#loading in data from other case study scripts
mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_acerifolia'
setwd(mydir)
load("combined_results_q_acerifolia.Rdata")

mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_engelmannii'
setwd(mydir)
load("combined_results_q_engelmannii.Rdata")

mydir = 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_oglethorpensis'
setwd(mydir)
load("combined_results_q_oglethorpensis.Rdata")

#preparing data by naming a column indicating species
q_acer = rep("Q. acerifolia", 200)
combined_q_acerifolia$species = q_acer

q_engel = rep("Q. engelmannii", 200)
combined_q_engelmannii$species = q_engel

q_ogle = rep("Q. oglethorpensis", 200)
combined_q_oglethorpensis$species = q_ogle

#combining all results into one larger dataframe for plotting on one graph
all_case_studies = rbind(combined_q_acerifolia, combined_q_engelmannii, combined_q_oglethorpensis)

#Creating graph with all case study species compared
p = ggplot(all_case_studies, aes(x=species, y=prop_all, fill=strategy)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ggtitle("Case study species") +
  xlab("Species") +
  ylab("Proportion of alleles captured") +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw() #+
  #theme(legend.position = "none")
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))
