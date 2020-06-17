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
for(i in 1:length(example_files)) {
  example_data[[i]] = arp2gen(example_files[i])
}
