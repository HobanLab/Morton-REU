setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\forest_extended")
library(diveRsity)
library(adegenet)

list.files(pattern = ".arp$")

#creating a list to store all file names
forest_files = list.files(pattern = ".arp$")
forest_files

#create an empty list to store data read in from files
forest_data = list()

#loop to read in data & create .gen files
for(i in 1:length(forest_files)) {
  forest_data[[i]] = arp2gen(forest_files[i])
}


#making genind object by reading in .gen files with loop
forest_genind = list()

#list of filenames
forest_files_gen = list.files(pattern = ".gen$")
forest_files_gen

#loop creates a list of separate genind objects
for(i in 1:length(forest_files_gen)){
  forest_genind[[i]] = read.genepop(forest_files_gen[i], ncode = 3)
}

#using repool function to combine list of genind objects to one 
forest_genind_complete = repool(forest_genind, list = FALSE)
forest_genind_complete

#converting from genind to genpop
forest_genpop = genind2genpop(forest_genind_complete)
