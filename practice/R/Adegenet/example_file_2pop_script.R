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

#not working
rep_1 = genind_list[1]
rep_1 #prints genind object
is.genind(rep_1)
rep_1@loc.n.all

#repool to one genind object
genind_complete = repool(genind_list, list = FALSE)
genind_complete

#summary
sum_genind = summary(genind_complete)

#par(mfrow=c(2,2))
plot(sum_genind$n.by.pop, sum_genind$pop.n.all, xlab="Sample size", 
     ylab="Number of alleles", main="Alleles numbers and sample sizes",
     type="n")
text(sum_genind$n.by.pop, sum_genind$pop.n.all, lab=names(sum_genind$n.by.pop))

barplot(sum_genind$loc.n.all, ylab="Number of alleles", main="Number of alleles per locus")
#heterozygosity
barplot(sum_genind$Hexp-sum_genind$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")

barplot(sum_genind$n.by.pop, main ="Sample sizes per population", ylab="Number of genotypes", las=3)



#number of alleles for each marker
genind_complete@loc.n.all
#number of loci
nLoc(genind_complete)
nAll(genind_complete)
#number of individuals in the object
nInd(genind_complete)
