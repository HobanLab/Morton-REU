setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file")
library(diveRsity)
library(adegenet)

list.files(pattern = ".arp$")

#creating a list to store all file names
example_files = list.files(pattern = ".arp$")
example_files

#create an empty list to store data read in from files
example_data = list()

#loop to read in data
#error here
#define new function to account for error
source("arp2gen_edited.R")

#loop to read in data to get .gen files
#runs but no .gen files created (could it have been placed somewhere else?)
for(i in 1:length(example_files)) {
  example_data[[i]] = arp2gen_edited(example_files[i])
}

##Converting to genind objects
example_file_genind = list()
example_filenames_gen = list.files(pattern = ".gen$")
example_filenames_gen

for(j in 1:length(example_filenames_gen)) {
  example_file_genind[[j]] = read.genepop(example_filenames_gen[j], ncode = 3)
}




