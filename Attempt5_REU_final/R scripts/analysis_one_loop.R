library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#root directory
#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt5_REU_final\\Simulations"
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

#import functions
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

##converting .arp to .gen for all combinations and scenarios
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    import_arp2gen_files(paste(my_dir,combinations[i],scenarios[j],sep=""), ".arp$")
  }
}

#****************************************************************************************************************************************************
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

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    for(k in 1:length(list_files)) {
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal[[j]],])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop[[j]],])>0)
      
      total_alleles = ncol(temp_genind@tab)
      
      if(i == 1) {
        results_highMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        total_alleles_highMig[j,k] = total_alleles
        
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
      } else {
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        total_alleles_lowMig[j,k] = total_alleles
        
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
      }
    }
  }
}