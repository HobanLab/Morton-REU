library(adegenet)
library(diveRsity)

setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop")

#import functions
import_arp2gen_files = function(mypath, mypattern) {
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

import_gen2genind_objects = function(mypath, mypattern) {
  temp_list_3 = list.files(mypath, mypattern)
  temp_list_4 = list(length = length(temp_list_3))
  for(j in 1:length(temp_list_3)){temp_list_4[[j]]=read.genepop(temp_list_3[j], ncode=3)}
  temp_list_4
}

import_arp2gen_files("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop", ".arp$")

genind_list = list()
genind_list = import_gen2genind_objects("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop", ".gen$")

rep_1 = genind_list[[1]]#note syntax
rep_1 #prints genind object
is.genind(rep_1)
rep_1@loc.n.all


#summary
sum_rep_1 = summary(rep_1)
sum_rep_1$pop.n.all

#plotting statistics for 
plot(sum_rep_1$n.by.pop, sum_rep_1$pop.n.all, xlab="Sample size", 
     ylab="Number of alleles", main="Alleles numbers and sample sizes")
text(sum_rep_1$n.by.pop, sum_rep_1$pop.n.all, lab=names(sum_rep_1$n.by.pop))

barplot(sum_rep_1$loc.n.all, ylab="Number of alleles", main="Number of alleles per locus")
#heterozygosity
barplot(sum_rep_1$Hexp-sum_rep_1$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")

barplot(sum_rep_1$n.by.pop, main ="Sample sizes per population", ylab="Number of genotypes", las=3)


