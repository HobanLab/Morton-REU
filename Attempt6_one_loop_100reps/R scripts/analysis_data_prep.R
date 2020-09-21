#analysis_data_prep
#this file converts the results arrays into data frames for utility in ggplot2
#columns are added to the dataframes to keep track of scenario and strategy variables

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt6_one_loop_100reps\\R scripts")
load("results_highMig_highSamp.Rdata")
load("results_lowMig_highSamp.Rdata")
load("results_highMig_lowSamp.Rdata")
load("results_lowMig_lowSamp.Rdata")

#PREPARING DATA (for plotting and statistical analyses)
#converting results arrays to dataframes
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
#variable to keep track of scenario
scenario = rep(c(1,2,3,4,5,6,7,8,9), 100)
#variables to keep track of strategy
equal_strategy = rep("equal", 900)
prop_strategy = rep("proportional", 900)

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

#combining dataframes (equal and proportional) for plotting
combined_highMig_highSamp = rbind(results_highMig_highSamp_equal_long, results_highMig_highSamp_prop_long)
combined_lowMig_highSamp = rbind(results_lowMig_highSamp_equal_long, results_lowMig_highSamp_prop_long)
combined_highMig_lowSamp = rbind(results_highMig_lowSamp_equal_long, results_highMig_lowSamp_prop_long)
combined_lowMig_lowSamp = rbind(results_lowMig_lowSamp_equal_long, results_lowMig_lowSamp_prop_long)

save(combined_highMig_highSamp, combined_lowMig_highSamp, combined_highMig_lowSamp, combined_lowMig_lowSamp, file="combined_dataframes.Rdata")

