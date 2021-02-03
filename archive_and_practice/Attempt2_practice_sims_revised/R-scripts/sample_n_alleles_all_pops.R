library(adegenet)
library(diveRsity)

my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations"

#list of scenarios
#simulation file folder directories
scenarios = c("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\equal_pop",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_extreme",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_moderate",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_strong",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_weak")

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
#5 populations
#10 replicates
#5 scenarios
results = array(0, dim = c(5,10,5))

#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(scenarios[i])
  list_files = list.files(path = scenarios[i], pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    #total number of alleles in the population
    total_n_alleles = ncol(temp_genind@tab)
    t = 10
    #number of alleles in 10 samples
    sample_n_alleles = sum(colSums(temp_genind@tab[sample(1:nrow(temp_genind@tab), t),])>0)
    results[,j,i] = sample_n_alleles
  }
}

#look at results
results
