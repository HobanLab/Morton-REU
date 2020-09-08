#analysis_conversions_calculations
#this script contains the conversions for files from Arlequin format to genepop format
#and loops that analyze the data involving converting genepop files to genind objects, 
#calculating the proportion of alleles captured, expected heterozygosity, and total number of alleles

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#defining root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt5_allcombinations_100reps\\Simulations"
setwd(my_dir)

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig", "\\lowMig")

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

#import function
#converts arlequin files to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

#loop converting .arp to .gen for all combinations and scenarios
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    import_arp2gen_files(paste(my_dir,combinations[i],scenarios[j],sep=""), ".arp$")
  }
}

#HIGH INTENSITY SAMPLING
#defining results arrays to store proportion of alleles captured
results_highMig_highSamp_equal = array(0, dim = c(9,100))
results_lowMig_highSamp_equal = array(0, dim = c(9,100))
results_highMig_highSamp_prop = array(0, dim = c(9,100))
results_lowMig_highSamp_prop = array(0, dim = c(9,100))

#defining arrays to store the total alleles captured
total_alleles_highMig = array(0, dim = c(9, 100))
total_alleles_lowMig = array(0, dim = c(9, 100))
#defining array to store the expected heterozygosity for each locus
hexp_highMig = array(0, dim = c(9,100))
hexp_lowMig = array(0, dim = c(9,100))

#creating list of vectors representing which individuals to sample from genind object
#sampling equally from each population - equal strategy
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

#sampling 10% from each population - proportional strategy
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

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  #loop going through each scenario for each combination
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #loop going through every replicate for each scenario
    for(k in 1:length(list_files)) {
      #creating a temporary genind object
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining how many alleles to sample based on rows_to_samp
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal[[j]],])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop[[j]],])>0)
      
      #saving the total alleles present
      total_alleles = ncol(temp_genind@tab)
      
      #if high migration 
      if(i == 1) {
        #saving proportion of alleles captured for both equal and proportional strategies
        results_highMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_highMig[j,k] = total_alleles
        
        #saving expected heterozygosity
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
      } else { #if low migration
        #saving proportion of alleles captured 
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_lowMig[j,k] = total_alleles
        
        #saving expected heterozyogisty
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
      }
    }
  }
}

#LOW INTENSITY SAMPLING
#defining results arrays to store proportion of alleles captured
results_highMig_lowSamp_equal = array(0, dim = c(9,100))
results_lowMig_lowSamp_equal = array(0, dim = c(9,100))
results_highMig_lowSamp_prop = array(0, dim = c(9,100))
results_lowMig_lowSamp_prop = array(0, dim = c(9,100))

#defining arrays to store the total alleles captured
total_alleles_highMig = array(0, dim = c(9, 100))
total_alleles_lowMig = array(0, dim = c(9, 100))
#defining array to store the expected heterozygosity for each locus
hexp_highMig = array(0, dim = c(9,100))
hexp_lowMig = array(0, dim = c(9,100))

#creating list of vectors representing which individuals to sample from genind object
#sampling equally from each population - equal strategy
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

#sampling 5% from each population
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

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    #going through every combinatio and every scenario
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #going through every replicate for each scenario
    for(k in 1:length(list_files)) {
      #creating a temporary genind object
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining a sample of alleles based on rows_to_samp
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal[[j]],])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop[[j]],])>0)
      
      #aclculating total alleles present
      total_alleles = ncol(temp_genind@tab)
      
      if(i == 1) {
        #saving proportion of alleles captured for both strategies
        results_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_highMig[j,k] = total_alleles
      } else {
        #saving proportion of alleles captured for both strategies
        results_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles 
        total_alleles_lowMig[j,k] = total_alleles
      }
    }
  }
}

#saving results
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt5_allcombinations_100reps\\R scripts")
save(results_highMig_highSamp_equal, results_highMig_highSamp_prop, file="results_highMig_highSamp.Rdata")
save(results_lowMig_highSamp_equal, results_lowMig_highSamp_prop, file="results_lowMig_highSamp.Rdata")
save(results_highMig_lowSamp_equal, results_highMig_lowSamp_prop, file="results_highMig_lowSamp.Rdata")
save(results_lowMig_lowSamp_equal, results_lowMig_lowSamp_prop, file="results_lowMig_lowSamp.Rdata")
