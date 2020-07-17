library(adegenet)
library(diveRsity)
library(ggplot2)
library(tidyr)

my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp"

#list of scenarios
#simulation file folder directories
scenarios = c("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen1_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen2_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen3_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen4_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen5_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen6_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen7_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen8_highMig_highSamp",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt3_full_factorial_lowRep\\Simulations\\highMig_highSamp\\scen9_highMig_highSamp")

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
  import_arp2gen_files(scenarios[i], ".arp$")
}

#creating results array to store the results
#9 scenarios (each with 5 populations)
#10 replicates
results_highMig_highSamp = array(0, dim = c(9,10))

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp = list(length = length(scenarios))
rows_to_samp[[1]] = c(sample(1:30,3), sample(31:130,10), sample(131:230,10), sample(231:330,10), sample(331:1500,117))
rows_to_samp[[2]] = c(sample(1:40,4), sample(41:190,15), sample(191:340,15), sample(341:490,15), sample(491:1500,101))
rows_to_samp[[3]] = c(sample(1:50,5), sample(51:250,20), sample(251:450,20), sample(451:650,20), sample(651:1500,85))
rows_to_samp[[4]] = c(sample(1:100,10), sample(101:300,20), sample(301:500,20), sample(501:700,20), sample(701:1500,80))
rows_to_samp[[5]] = c(sample(1:150,15), sample(151:350,20), sample(351:550,20), sample(551:750,20), sample(750:1500,75))
rows_to_samp[[6]] = c(sample(1:200,20), sample(201:450,25), sample(451:700,25), sample(701:950,25), sample(951:1500,55))
rows_to_samp[[7]] = c(sample(1:200,20), sample(201:500,30), sample(501:800,30), sample(801:1100,30), sample(1101:1500,40))
rows_to_samp[[8]] = c(sample(1:290,29), sample(291:590,30), sample(591:890,30), sample(891:1190,30), sample(1191:1500,31))
rows_to_samp[[9]] = c(sample(1:300,30), sample(301:600,30), sample(601:900,30), sample(901:1200,30), sample(1201:1500,30))


#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(scenarios[i])
  list_files = list.files(path = scenarios[i], pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp[[i]],])>0)
    total_alleles = ncol(temp_genind@tab)
    results_highMig_highSamp[i,j] = sample_n_alleles/total_alleles
  }
}

#look at results
round(results_highMig_highSamp, 3)

#plotting
#this is a mess
results_plot = as.data.frame(results_highMig_highSamp) #converting to a dataframe
results_plot_long = gather(results_plot, replicate, prop_all) #converting to 'long' format
scenario = rep(c(1,2,3,4,5,6,7,8,9), 10) #repeating the number of scenarios x number replicates (10)
factor(scenario) #making this a factor for plotting purposes
results_plot_long$scenario=scenario #creating column to keep track of scenarios

ggplot(results_plot_long, aes(x=scenario, y=prop_all, group=scenario, fill=factor(scenario))) + 
  geom_boxplot() +
  ggtitle("Alleles captured across scenarios: equal strategy") +
  ylim(0.85,1) +
  xlim(0,10)
