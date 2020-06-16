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
list.filenames

#create an empty list to store data read in from files
list.data = list()

#loop to read in data
for(i in 1:length(list.filenames)){
  list.data[[i]] = arp2gen(list.filenames[i])
}

#adding names to data
names(list.data) = list.filenames

list.data

