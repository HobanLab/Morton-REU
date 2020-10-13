library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(tidyr)

#root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt4_full_factorial_100rep\\Simulations"
setwd(my_dir)

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig", "\\lowMig")

#list of scenarios
#simulation sub-folder directories
scenarios_highMig = c("\\scen1_highMig_highSamp",
                   "\\scen2_highMig_highSamp",
                   "\\scen3_highMig_highSamp",
                   "\\scen4_highMig_highSamp",
                   "\\scen5_highMig_highSamp",
                   "\\scen6_highMig_highSamp",
                   "\\scen7_highMig_highSamp",
                   "\\scen8_highMig_highSamp",
                   "\\scen9_highMig_highSamp")
scenarios_lowMig = c("\\scen1_lowMig_highSamp",
                   "\\scen2_lowMig_highSamp",
                   "\\scen3_lowMig_highSamp",
                   "\\scen4_lowMig_highSamp",
                   "\\scen5_lowMig_highSamp",
                   "\\scen6_lowMig_highSamp",
                   "\\scen7_lowMig_highSamp",
                   "\\scen8_lowMig_highSamp",
                   "\\scen9_lowMig_highSamp")

setwd(paste(my_dir, combinations[1]))
