#analysis_graphics_results
#this script contains code for creating plots, and Wilcoxon ranks sums tests
#we saved the p-values from the tests in a matrix

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt6_complete_loop\\R scripts")
load("results_highMig_highSamp.Rdata")
load("results_lowMig_highSamp.Rdata")
load("results_highMig_lowSamp.Rdata")
load("results_lowMig_lowSamp.Rdata")
load("combined_dataframes.Rdata")

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig_highSamp", "\\lowMig_highSamp", "\\highMig_lowsamp", "\\lowMig_lowSamp")

#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1",
              "\\scen2",
              "\\scen3",
              "\\scen4",
              "\\scen5",
              "\\scen6",
              "\\scen7",
              "\\scen8",
              "\\scen9")

#PLOTTING RESULTS
#High migration high sampling
p = ggplot(combined_highMig_highSamp, aes(x=factor(scenario), y=prop_all, fill=strategy, color=factor(scenario))) + 
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE, label.y = c(0.96,0.97,0.98,0.98,0.985)) +
  ggtitle("High migration high sampling intensity") +
  xlab("Scenarios") +
  ylab("Proportion of alleles captured") +
  labs(color = "Scenario", fill = "Sample strategy") +
  ylim(0.85,1) +
  scale_fill_brewer() +
  theme_bw()
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))

#Low migration high sampling
p = ggplot(combined_lowMig_highSamp, aes(x=factor(scenario), y=prop_all, fill=strategy, color=factor(scenario))) + 
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE, label.y = c(0.985,0.985,0.985,0.985,0.99, 0.99)) +
  ggtitle("Low migration high sampling intensity") +
  xlab("Scenarios") +
  ylab("Proportion of alleles captured") +
  labs(color = "Scenario", fill = "Sample strategy") +
  ylim(0.85,1) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw()
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))

#High migration low sampling
p = ggplot(combined_highMig_lowSamp, aes(x=factor(scenario), y=prop_all, fill=strategy, color=factor(scenario))) + 
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE, label.y = c(0.94,0.94,0.95,0.96,0.96)) +
  ggtitle("High migration low sampling intensity") +
  xlab("Scenarios") +
  ylab("Proportion of alleles captured") +
  labs(color = "Scenario", fill = "Sample strategy") +
  ylim(0.85,1) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw()
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))

#Low migration low sampling
p = ggplot(combined_lowMig_lowSamp, aes(x=factor(scenario), y=prop_all, fill=strategy, color=factor(scenario))) + 
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE, label.y = c(0.96,0.97,0.97,0.975,0.975,0.975)) +
  ggtitle("Low migration low sampling intensity") +
  xlab("Scenarios") +
  ylab("Proportion of alleles captured") +
  labs(color = "Scenario", fill = "Sample strategy") +
  ylim(0.85,1) +
  scale_fill_brewer() +
  theme(axis.text=element_text(size=30, face="bold")) +
  theme_bw()
p + theme(axis.text = element_text(size = 11, face = "bold"), axis.title = element_text(size = 14))

#STATISTICAL ANALYSES
#Wilcoxon rank sums test
#vector to store the p-values of each test for each combination
p_val_highMig_highSamp = c(rep(0, 9))
p_val_lowMig_highSamp = c(rep(0, 9))
p_val_highMig_lowSamp = c(rep(0, 9))
p_val_lowMig_lowSamp = c(rep(0, 9))

#loop to calculate Wilcox tests and save p-values in vectors
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    if(i == 1) { #high mig high samp
      x_var = results_highMig_highSamp_prop[j,]
      y_var = results_highMig_highSamp_equal[j,]
      test_result = wilcox.test(x_var, y_var)
      p_val_highMig_highSamp[j] = test_result$p.value
    } else if(i == 2) { #low mig high samp
      x_var = results_lowMig_highSamp_prop[j,]
      y_var = results_lowMig_highSamp_equal[j,]
      test_result = wilcox.test(x_var, y_var)
      p_val_lowMig_highSamp[j] = test_result$p.value
    } else if(i == 3) { #high mig low samp
      x_var = results_highMig_lowSamp_prop[j,]
      y_var = results_highMig_lowSamp_equal[j,]
      test_result = wilcox.test(x_var, y_var)
      p_val_highMig_lowSamp[j] = test_result$p.value
    } else { #low mig low samp
      x_var = results_lowMig_lowSamp_prop[j,]
      y_var = results_lowMig_lowSamp_equal[j,]
      test_result = wilcox.test(x_var, y_var)
      p_val_lowMig_lowSamp[j] = test_result$p.value
    }
  }
}

#merging all vectors together into a matrix combining all p values across combinations
p_value_matrix = matrix(0, nrow = 4, ncol = 9)
p_value_matrix[1,] = p_val_highMig_highSamp
p_value_matrix[2,] = p_val_lowMig_highSamp
p_value_matrix[3,] = p_val_highMig_lowSamp
p_value_matrix[4,] = p_val_lowMig_lowSamp

#naming rows and columns for matrix and printing results
rownames(p_value_matrix) = c("highMig_highSamp", "lowMig_highSamp", "highMig_lowSamp", "lowMig_lowSamp")
colnames(p_value_matrix) = c("scenario 1", "scenario 2", "scenario 3", "scenario 4", "scenario 5", "scenario 6", "scenario 7", "scenario 8", "scenario 9")
matrix(round(p_value_matrix, 5), nrow=4, ncol=9)

#p adjustment methods
matrix(round(p.adjust(p_value_matrix, method = "BH"), 5), nrow=4, ncol=9)
matrix(round(p.adjust(p_value_matrix, method = "BY"), 5), nrow=4, ncol=9)
matrix(round(p.adjust(p_value_matrix, method = "bonferroni"), 5), nrow=4, ncol=9)
