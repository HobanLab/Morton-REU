#Script to analyze 2 different scenarios
#scenario 1: all equal size populations
#scenario 2: unequal size populations
#script imports data, analyzes, and saves results to an array

library(adegenet)
library(diveRsity)

my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt_1\\Simulations"

scenarios = c("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt_1\\Simulations\\all_same_size", 
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt_1\\Simulations\\all_diff_size")

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
import_arp2gen_files("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt_1\\Simulations\\all_same_size", ".arp$")
import_arp2gen_files("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt_1\\Simulations\\all_diff_size", ".arp$")

#creating results array to store the results
#5 populations
#10 replicates
#2 scenarios
results = array(0, dim = c(5,10,2))

#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(scenarios[i])
  list_files = list.files(path = scenarios[i], pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    #run summary
    temp_summary = summary(temp_genind)
    #extract statistics & store results
    results[,j,i] = temp_summary$pop.n.all
  }
}

#look at results
results


#Testing space
temp_genind = read.genepop(list_files[[1]], ncode=3)
temp_summary = summary(temp_genind)
temp_summary
temp_summary$pop.n.all
