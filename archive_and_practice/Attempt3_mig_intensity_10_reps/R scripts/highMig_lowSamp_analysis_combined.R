library(adegenet)
library(diveRsity)
library(ggplot2)
library(tidyr)

#root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_10rep\\Simulations\\highMig_lowSamp"
setwd(my_dir)


#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1_highMig_lowSamp",
              "\\scen2_highMig_lowSamp",
              "\\scen3_highMig_lowSamp",
              "\\scen4_highMig_lowSamp",
              "\\scen5_highMig_lowSamp",
              "\\scen6_highMig_lowSamp",
              "\\scen7_highMig_lowSamp",
              "\\scen8_highMig_lowSamp",
              "\\scen9_highMig_lowSamp")

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
results_highMig_lowSamp_equal = array(0, dim = c(9,10))

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp_equal = list(length = length(scenarios))
rows_to_samp_equal[[1]] = c(sample(1:30,15), sample(31:130,15), sample(131:230,15), sample(231:330,15), sample(331:1500,15))
rows_to_samp_equal[[2]] = c(sample(1:40,15), sample(41:190,15), sample(191:340,15), sample(341:490,15), sample(491:1500,15))
rows_to_samp_equal[[3]] = c(sample(1:50,15), sample(51:250,15), sample(251:450,15), sample(451:650,15), sample(651:1500,15))
rows_to_samp_equal[[4]] = c(sample(1:100,15), sample(101:300,15), sample(301:500,15), sample(501:700,15), sample(701:1500,15))
rows_to_samp_equal[[5]] = c(sample(1:150,15), sample(151:350,15), sample(351:550,15), sample(551:750,15), sample(750:1500,15))
rows_to_samp_equal[[6]] = c(sample(1:200,15), sample(201:450,15), sample(451:700,15), sample(701:950,15), sample(951:1500,15))
rows_to_samp_equal[[7]] = c(sample(1:200,15), sample(201:500,15), sample(501:800,15), sample(801:1100,15), sample(1101:1500,15))
rows_to_samp_equal[[8]] = c(sample(1:290,15), sample(291:590,15), sample(591:890,15), sample(891:1190,15), sample(1191:1500,15))
rows_to_samp_equal[[9]] = c(sample(1:300,15), sample(301:600,15), sample(601:900,15), sample(901:1200,15), sample(1201:1500,15))


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
    results_highMig_lowSamp_equal[i,j] = sample_n_alleles/total_alleles
  }
}

#look at results
round(results_highMig_lowSamp_equal, 3)

#**************************************************************************************************************************************************************
#proportional strategy
#creating results array to store the results for the proportional strategy
#9 scenarios (each with 5 populations)
#10 replicates
results_highMig_lowSamp_prop = array(0, dim = c(9,10))

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp_prop = list(length = length(scenarios))
rows_to_samp_prop[[1]] = c(sample(1:30,2), sample(31:130,5), sample(131:230,5), sample(231:330,5), sample(331:1500,59))
rows_to_samp_prop[[2]] = c(sample(1:40,2), sample(41:190,8), sample(191:340,8), sample(341:490,8), sample(491:1500,51))
rows_to_samp_prop[[3]] = c(sample(1:50,3), sample(51:250,10), sample(251:450,10), sample(451:650,10), sample(651:1500,43))
rows_to_samp_prop[[4]] = c(sample(1:100,5), sample(101:300,10), sample(301:500,10), sample(501:700,10), sample(701:1500,40))
rows_to_samp_prop[[5]] = c(sample(1:150,8), sample(151:350,10), sample(351:550,10), sample(551:750,10), sample(750:1500,38))
rows_to_samp_prop[[6]] = c(sample(1:200,10), sample(201:450,13), sample(451:700,13), sample(701:950,13), sample(951:1500,28))
rows_to_samp_prop[[7]] = c(sample(1:200,10), sample(201:500,15), sample(501:800,15), sample(801:1100,15), sample(1101:1500,20))
rows_to_samp_prop[[8]] = c(sample(1:290,15), sample(291:590,15), sample(591:890,15), sample(891:1190,15), sample(1191:1500,15))
rows_to_samp_prop[[9]] = c(sample(1:300,15), sample(301:600,15), sample(601:900,15), sample(901:1200,15), sample(1201:1500,15))

#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(paste(my_dir,scenarios[i],sep=""))
  list_files = list.files(path = paste(my_dir,scenarios[i],sep=""), pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    #check population orer
    print(table(temp_genind@pop))
    #sampling alelles from the population
    sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp_prop[[i]],])>0)
    #keeping track of the total alleles
    total_alleles = ncol(temp_genind@tab)
    #calculating proportion and saving the results
    results_highMig_lowSamp_prop[i,j] = sample_n_alleles/total_alleles
  }
}

#look at results
round(results_highMig_lowSamp_prop, 3)

#***********************************************************************************************************************************************************
#converting results arrays to matrices
#Equal results array conversion:
results_plot_equal = as.data.frame(results_highMig_lowSamp_equal)
results_plot_equal_long = gather(results_plot_equal, replicate, prop_all)
scenario = rep(c(1,2,3,4,5,6,7,8,9), 10) #repeating the number of scenarios x number replicates (10)
factor(scenario)
results_plot_equal_long$scenario = scenario
strategy = rep("equal", 90)#repeat equal 90 times to keep track that this is the equal strategy
results_plot_equal_long$strategy=strategy

#Proportional results array conversion:
results_plot_prop = as.data.frame(results_highMig_lowSamp_prop)
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
  ggtitle("High migration/low sampling") +
  ylim(0.85,1) +
  scale_fill_brewer()
#note:can take out color = to get rid of mess

#faceted results
#separate plots for each scenario
ggplot(combined_results, aes(x=factor(scenario), y=prop_all, fill=strategy)) + 
  geom_boxplot() +
  ggtitle("High migration/low sampling") +
  ylim(0.85,1) +
  facet_wrap(~scenario, scale="free") +
  scale_fill_brewer(palette = "blues")
