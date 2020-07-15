library(adegenet)
library(diveRsity)
library(ggplot2)
library(tidyr)

#root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp"
setwd(my_dir)


#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1_highMig_highSamp",
              "\\scen2_highMig_highSamp",
              "\\scen3_highMig_highSamp",
              "\\scen4_highMig_highSamp",
              "\\scen5_highMig_highSamp",
              "\\scen6_highMig_highSamp",
              "\\scen7_highMig_highSamp",
              "\\scen8_highMig_highSamp",
              "\\scen9_highMig_highSamp")

#import functions
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

import_gen2genind_objects = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_3 = list.files(mypath, mypattern)
  temp_list_4 = list(length = length(temp_list_3))
  for(j in 1:length(temp_list_3)){temp_list_4[[j]]=read.genepop(temp_list_3[j], ncode=3)}
  temp_list_4
}

##converting .arp to .gen
for(i in 1:length(scenarios)) {
  import_arp2gen_files(paste(my_dir,scenarios[i],sep=""), ".arp$")
}

#****************************************************************************************************************************************************
#equal strategy
#creating results array to store the results for the equal strategy
#9 scenarios (each with 5 populations)
#10 replicates
results_highMig_highSamp_equal = array(0, dim = c(9,10))

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp_equal = list(length = length(scenarios))
rows_to_samp_equal[[1]] = c(sample(1:30,30), sample(31:130,30), sample(131:230,30), sample(231:330,30), sample(331:1500,30))
rows_to_samp_equal[[2]] = c(sample(1:40,30), sample(41:190,30), sample(191:340,30), sample(341:490,30), sample(491:1500,30))
rows_to_samp_equal[[3]] = c(sample(1:50,30), sample(51:250,30), sample(251:450,30), sample(451:650,30), sample(651:1500,30))
rows_to_samp_equal[[4]] = c(sample(1:100,30), sample(101:300,30), sample(301:500,30), sample(501:700,30), sample(701:1500,30))
rows_to_samp_equal[[5]] = c(sample(1:150,30), sample(151:350,30), sample(351:550,30), sample(551:750,30), sample(750:1500,30))
rows_to_samp_equal[[6]] = c(sample(1:200,30), sample(201:450,30), sample(451:700,30), sample(701:950,30), sample(951:1500,30))
rows_to_samp_equal[[7]] = c(sample(1:200,30), sample(201:500,30), sample(501:800,30), sample(801:1100,30), sample(1101:1500,30))
rows_to_samp_equal[[8]] = c(sample(1:290,30), sample(291:590,30), sample(591:890,30), sample(891:1190,30), sample(1191:1500,30))
rows_to_samp_equal[[9]] = c(sample(1:300,30), sample(301:600,30), sample(601:900,30), sample(901:1200,30), sample(1201:1500,30))


#loop through scenarios
#equal strategy
for(i in 1:length(scenarios)) {
  setwd(paste(my_dir,scenarios[i],sep=""))
  list_files = list.files(path = paste(my_dir,scenarios[i],sep=""), pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    #sampling alleles from each population
    sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp_equal[[i]],])>0)
    #total alleles
    total_alleles = ncol(temp_genind@tab)
    #saving results
    results_highMig_highSamp_equal[i,j] = sample_n_alleles/total_alleles
  }
}

#look at results
round(results_highMig_highSamp_equal, 3)

#**************************************************************************************************************************************************************
#proportional strategy
#creating results array to store the results for the proportional strategy
#9 scenarios (each with 5 populations)
#10 replicates
results_highMig_highSamp_prop = array(0, dim = c(9,10))

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp_prop = list(length = length(scenarios))
rows_to_samp_prop[[1]] = c(sample(1:30,3), sample(31:130,10), sample(131:230,10), sample(231:330,10), sample(331:1500,117))
rows_to_samp_prop[[2]] = c(sample(1:40,4), sample(41:190,15), sample(191:340,15), sample(341:490,15), sample(491:1500,101))
rows_to_samp_prop[[3]] = c(sample(1:50,5), sample(51:250,20), sample(251:450,20), sample(451:650,20), sample(651:1500,85))
rows_to_samp_prop[[4]] = c(sample(1:100,10), sample(101:300,20), sample(301:500,20), sample(501:700,20), sample(701:1500,80))
rows_to_samp_prop[[5]] = c(sample(1:150,15), sample(151:350,20), sample(351:550,20), sample(551:750,20), sample(750:1500,75))
rows_to_samp_prop[[6]] = c(sample(1:200,20), sample(201:450,25), sample(451:700,25), sample(701:950,25), sample(951:1500,55))
rows_to_samp_prop[[7]] = c(sample(1:200,20), sample(201:500,30), sample(501:800,30), sample(801:1100,30), sample(1101:1500,40))
rows_to_samp_prop[[8]] = c(sample(1:290,29), sample(291:590,30), sample(591:890,30), sample(891:1190,30), sample(1191:1500,31))
rows_to_samp_prop[[9]] = c(sample(1:300,30), sample(301:600,30), sample(601:900,30), sample(901:1200,30), sample(1201:1500,30))

#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(scenarios[i])
  list_files = list.files(path = scenarios[i], pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    #sampling alelles from the population
    sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp_prop[[i]],])>0)
    #keeping track of the total alleles
    total_alleles = ncol(temp_genind@tab)
    #calculating proportion and saving the results
    results_highMig_highSamp_prop[i,j] = sample_n_alleles/total_alleles
  }
}

#look at results
round(results_highMig_highSamp_prop, 3)

#***********************************************************************************************************************************************************
#converting results arrays to matrices
#Equal results array conversion:
results_plot_equal = as.data.frame(results_highMig_highSamp_equal)
results_plot_equal_long = gather(results_plot_equal, replicate, prop_all)
scenario = rep(c(1,2,3,4,5,6,7,8,9), 10) #repeating the number of scenarios x number replicates (10)
factor(scenario)
results_plot_equal_long$scenario = scenario
strategy = rep("equal", 90)#repeat equal 90 times to keep track that this is the equal strategy
results_plot_equal_long$strategy=strategy

#Proportional results array conversion:
results_plot_prop = as.data.frame(results_highMig_highSamp_prop)
results_plot_prop_long = gather(results_plot_prop, replicate, prop_all)
results_plot_prop_long$scenario = scenario
strategy = rep("proportional", 90)#repeat proportional to keep track that this is proportional strategy
results_plot_prop_long$strategy = strategy

#************************************************************************************************************************************************************
#merging the two data frames
combined_results = rbind(results_plot_equal_long, results_plot_prop_long)

#************************************************************************************************************************************************************
#plotting combined results
#all results on one plot
ggplot(combined_results, aes(x=factor(scenario), y=prop_all, fill=strategy, color=factor(scenario))) + 
  geom_boxplot() +
  ggtitle("High migration/high sampling") +
  ylim(0.85,1) +
  scale_fill_brewer()
#note:can take out color = to get rid of mess

#faceted results
#separate plots for each scenario
ggplot(combined_results, aes(x=factor(scenario), y=prop_all, fill=strategy)) + 
  geom_boxplot() +
  ggtitle("High migration/high sampling") +
  ylim(0.85,1) +
  facet_wrap(~scenario, scale="free") +
  scale_fill_brewer(palette = "blues")
 