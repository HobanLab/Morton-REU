library(adegenet)
library(diveRsity)

my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations"

#list of scenarios
#simulation file folder directories
scenarios = c("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\equal_pops",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_extreme",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_strong",
              "C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt2_fullpop\\Simulations\\unequal_pop_moderate",
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

#creating list of vectors representing rows to sample from genind object
#sampling 10% from each population
#5 pops in each scenario
rows_to_samp = list(length = length(scenarios))
rows_to_samp[[1]] = c(sample(1:300,30), sample(301:600,30), sample(601:900,30), sample(901:1200,30), sample(1201:1500,30))
rows_to_samp[[2]] = c(sample(1:20,2), sample(21:120,10), sample(121:220,10), sample(221:320, 10), sample(321:1500,118))
rows_to_samp[[3]] = c(sample(1:50,5), sample(51:250,20), sample(251:450,20), sample(451:650,20), sample(651:1500,85))
rows_to_samp[[4]] = c(sample(1:150,15), sample(151:350,20), sample(351:550,20), sample(551:750,20), sample(751:1500,75))
rows_to_samp[[5]] = c(sample(1:290,29), sample(291:590,30), sample(591:890,30), sample(891:1190,30), sample(1191:1500,31))


#loop through scenarios
for(i in 1:length(scenarios)) {
  setwd(scenarios[i])
  list_files = list.files(path = scenarios[i], pattern = ".gen$")
  #loop through replicates
  for(j in 1:length(list_files)) {
    #convert to genind
    temp_genind = read.genepop(list_files[[j]], ncode=3)
    sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp[[i]],])>0)
    results[,j,i] = sample_n_alleles
  }
}

#look at results
results

barplot(results[,,1], main="Equal populations", xlab="Replicates", ylab="sample_n_alleles")
barplot(results[,,2], main="Unequal populations: extreme", xlab="Replicates", ylab="sample_n_alleles")
barplot(results[,,3], main="Unequal populations: strong", xlab="Replicates", ylab="sample_n_alleles")
barplot(results[,,4], main="Unequal populations: moderate", xlab="Replicates", ylab="sample_n_alleles")
barplot(results[,,5], main="Unequal populations: weak", xlab="Replicates", ylab="sample_n_alleles")
