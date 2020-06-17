setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\forest_extended")
library(diveRsity)
library(adegenet)

list.files(pattern = ".arp$")

#creating a list to store all file names
forest_files = list.files(pattern = ".arp$")
forest_files

#create an empty list to store data read in from files
forest_data = list()

#loop to read in data
for(i in 1:length(forest_files)) {
  forest_data[[i]] = arp2gen(forest_files[i])
}
