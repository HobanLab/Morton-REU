#analysis_data_prep
#Code written by XXX XXX, XXX XXX, and XXX XXX XXX in collaboration
#This file converts the results arrays into data frames for utility in ggplot2
#Columns are added to the dataframes to keep track of scenario and strategy variables

#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#set working directory and load in data from calculations R script
setwd("C:\\Users\\kayle\\Documents\\XXX-XXX\\samp_pop_sims\\R scripts")
load("results_highMig_highSamp.Rdata")
load("results_lowMig_highSamp.Rdata")
load("results_highMig_lowSamp.Rdata")
load("results_lowMig_lowSamp.Rdata")

#PREPARING DATA (for plotting and statistical analyses)
#converting results arrays to dataframes (dataframes are easier to use in ggplot2 than matrices)
#we have 8 matrices total, holding results for which strategy (equal or prop.), migration rate, and sampling intensity
#***data frames can easily be visualized in R with the View() function
results_highMig_highSamp_equal_df = as.data.frame(results_highMig_highSamp_equal)
results_highMig_highSamp_prop_df = as.data.frame(results_highMig_highSamp_prop)
results_lowMig_highSamp_equal_df = as.data.frame(results_lowMig_highSamp_equal)
results_lowMig_highSamp_prop_df = as.data.frame(results_lowMig_highSamp_prop)
results_highMig_lowSamp_equal_df = as.data.frame(results_highMig_lowSamp_equal)
results_highMig_lowSamp_prop_df = as.data.frame(results_highMig_lowSamp_prop)
results_lowMig_lowSamp_equal_df = as.data.frame(results_lowMig_lowSamp_equal)
results_lowMig_lowSamp_prop_df = as.data.frame(results_lowMig_lowSamp_prop)

#converting results to 'long' format using gather()
results_highMig_highSamp_equal_long = gather(results_highMig_highSamp_equal_df, replicate, prop_all)
results_highMig_highSamp_prop_long = gather(results_highMig_highSamp_prop_df, replicate, prop_all)
results_lowMig_highSamp_equal_long = gather(results_lowMig_highSamp_equal_df, replicate, prop_all)
results_lowMig_highSamp_prop_long = gather(results_lowMig_highSamp_prop_df, replicate, prop_all)
results_highMig_lowSamp_equal_long = gather(results_highMig_lowSamp_equal_df, replicate, prop_all)
results_highMig_lowSamp_prop_long = gather(results_highMig_lowSamp_prop_df, replicate, prop_all)
results_lowMig_lowSamp_equal_long = gather(results_lowMig_lowSamp_equal_df, replicate, prop_all)
results_lowMig_lowSamp_prop_long = gather(results_lowMig_lowSamp_prop_df, replicate, prop_all)

#creating new columns for scenario and strategy
#variable to keep track of scenario - repeating 100 times for each replicate
scenario = rep(c(1,2,3,4,5,6,7,8,9), 100)
#variables to keep track of strategy
#these are each repeated 900 times since for each strategy, there is 9 scenarios and 100 replicates of each (900 total instances)
equal_strategy = rep("equal", 900)
prop_strategy = rep("proportional", 900)

#naming columns on existing dataframes to keep track of these variables when plotting
results_highMig_highSamp_equal_long$scenario = scenario
results_highMig_highSamp_prop_long$scenario = scenario
results_highMig_highSamp_equal_long$strategy = equal_strategy
results_highMig_highSamp_prop_long$strategy = prop_strategy

results_lowMig_highSamp_equal_long$scenario = scenario
results_lowMig_highSamp_prop_long$scenario = scenario
results_lowMig_highSamp_equal_long$strategy = equal_strategy
results_lowMig_highSamp_prop_long$strategy = prop_strategy

results_highMig_lowSamp_equal_long$scenario = scenario
results_highMig_lowSamp_prop_long$scenario = scenario
results_highMig_lowSamp_equal_long$strategy = equal_strategy
results_highMig_lowSamp_prop_long$strategy = prop_strategy

results_lowMig_lowSamp_equal_long$scenario = scenario
results_lowMig_lowSamp_prop_long$scenario = scenario
results_lowMig_lowSamp_equal_long$strategy = equal_strategy
results_lowMig_lowSamp_prop_long$strategy = prop_strategy

#combining dataframes vertically using rbind() (equal and proportional) for plotting
#then, both strategies can be visualized on the same plot
combined_highMig_highSamp = rbind(results_highMig_highSamp_equal_long, results_highMig_highSamp_prop_long)
combined_lowMig_highSamp = rbind(results_lowMig_highSamp_equal_long, results_lowMig_highSamp_prop_long)
combined_highMig_lowSamp = rbind(results_highMig_lowSamp_equal_long, results_highMig_lowSamp_prop_long)
combined_lowMig_lowSamp = rbind(results_lowMig_lowSamp_equal_long, results_lowMig_lowSamp_prop_long)

#saving data as Rdata file
save(combined_highMig_highSamp, combined_lowMig_highSamp, combined_highMig_lowSamp, combined_lowMig_lowSamp, file="combined_dataframes.Rdata")

