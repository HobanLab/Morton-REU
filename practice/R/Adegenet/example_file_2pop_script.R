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
pop_1 = genind_list[1]
pop_1
pop_1@loc.n.all

#repool to one genind object
genind_complete = repool(genind_list, list = FALSE)
genind_complete

#summary
summary(genind_complete)


#number of alleles for each marker
genind_complete@loc.n.all
#number of loci
nLoc(genind_complete)
nAll(genind_complete)
#number of individuals in the object
nInd(genind_complete)
