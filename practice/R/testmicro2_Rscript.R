setwd('C:\\Users\\kayle\\Documents\\Morton REU\\Morton-REU-master\\Simulations\\simcoal\\testmicro2')
library(diveRsity)
library(adegenet)

##Basic steps:
#Simulation -> Arlequin output files (.arp) -> convert Arlequin output file to genepop (.gen) format 
#-> load genepop format into R object genind -> do analysis

#lists all files in current directory with .arp extension
list.files(pattern = ".arp$")

#creating a list to store all file names
list.filenames = list.files(pattern = ".arp$")

#create an empty list to store data read in from files
list.data = list()

#loop to read in data
for(i in 1:length(list.filenames)){
  list.data[[i]] = arp2gen(list.filenames[i])
}

#adding names to data
#names(list.data) = list.filenames
#list.data
#renaming - genepop format
testmicro2_genepop = list.data

#*********************************************
#genepop to genind
#another loop to convert?
#empty list to store
testmicro2_genind = list()

list.files(pattern = ".gen$")
list.filenames.gen = list.files(pattern = ".gen$")

#loop to read in .gen files to form genind object
for(i in 1:length(list.filenames.gen)){
  testmicro2_genind[[i]] = read.genepop(list.filenames.gen[i], ncode = 3)
}


#attempting to make a function out of it
#import.arp2gen.files = function(mypath, mypattern, ...) {
  #tmp.list.1 = list.files(mypath, pattern = mypattern)
  #tmp.list.2 = list(length = length(tmp.list.1))
  #for(i in 1:length(tmp.list.1)){tmp.list.2[[i]]=arp2gen(tmp.list.1[i], ...)}
  #names(tmp.list.2) = tmp.list.1
  #tmp.list.2
#}
#testmicro2_import = import.arp2gen.files("C:\\Users\\kayle\\Documents\\Morton REU\\Morton-REU-master\\Simulations\\simcoal\\testmicro2", ".arp$")
