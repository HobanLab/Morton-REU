library(adegenet)
library(diveRsity)

setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\all_EU_lmlb")

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

import_arp2gen_files("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\all_EU_lmlb", ".arp")

all_EU_lmlb_genind = import_gen2genind_objects("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\all_EU_lmlb", ".gen$")

all_EU_lmlb_genind[[1]]
