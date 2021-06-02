#Code written by Emily Schumacher and Dr. Sean Hoban in collaboration

####################################
########## Load Libraries ##########
####################################

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

#####################################
###### Load Directories #############
#####################################

#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims_bottlenecks\\bottleneck_2\\Simulations"
setwd(my_dir)

###allelic capture functions 
source("C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims_bottlenecks\\bottleneck_2\\R-scripts\\Fa_sample_funcs.R")
##functions
colMax <- function(data) sapply(data, max, na.rm = TRUE)
sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
  p<-length(vect_pop_ID)
  if (p>1) {
    for (p in 1:length(vect_pop_ID))
      alleles[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)
    alleles<-colSums(alleles)
  } else {alleles<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)}
  
  alleles
}    

########################
######## Lists #########
########################

##allele frequency categories 
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#list of combinations
#combination sub-folder directories
#here, we are using the files for 'high sampling intensity'
#but this doesn't make a difference. low sampling intensity files
#are identical copies, just copied for easier naviagtion
#since they are the same, we are using the data from _hgihSamp
combinations = c("\\highMig_highSamp", "\\lowMig_highSamp")

#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1",
              "\\scen2",
              "\\scen3",
              "\\scen4",
              "\\scen5",
              "\\scen6",
              "\\scen7",
              "\\scen8",
              "\\scen9")

#################################################################
##### Run Sampling Code and Determine Allelic Capture ###########
#################################################################

##matrices
highmig_all_existing_by_sp_reps <- array(dim = c(100,9,9))
lowmig_all_existing_by_sp_reps <- array(dim = c(100,9,9))

##intermediate
highmig_equal_alleles <- array(dim = c(100,9,9))
highmig_prop_alleles <- array(dim = c(100,9,9))

highmig_equal_per <- array(dim = c(100,9,9))
highmig_prop_per <- array(dim = c(100,9,9))

##outputs
highmig_alleles_existing_by_cat <- matrix(nrow = 9, ncol = 9)
highmig_equal_mean_all_cap <- matrix(nrow = 9, ncol = 9)
highmig_prop_mean_all_cap <- matrix(nrow = 9, ncol = 9)
highmig_equal_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)
highmig_prop_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)

##intermediate
lowmig_equal_alleles <- array(dim = c(100,9,9))
lowmig_prop_alleles <- array(dim = c(100,9,9))

lowmig_equal_per <- array(dim = c(100,9,9))
lowmig_prop_per <- array(dim = c(100,9,9))

##outputs
lowmig_alleles_existing_by_cat <- matrix(nrow = 9, ncol = 9)
lowmig_equal_mean_all_cap <- matrix(nrow = 9, ncol = 9)
lowmig_prop_mean_all_cap <- matrix(nrow = 9, ncol = 9)
lowmig_equal_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)
lowmig_prop_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)

#This script will also collect in the freq_v_cap lists: a list of all the alleles in every simulation, their frequency in the population, and how many times they are captured (which also identifies those captured 0 times)

#Lists to hold these matrices of alleles and capture counts
#Each list element is a scenario, 1 to 9
#inside each list element is a matrix, with each row being an individual allele 
#col 1 = number of counts of that allele captured, col 2 = that alleles frequency in the simulated data (in situ)
#alleles from every replicate will all be concatenated with rbind into a super long list 
freq_v_cap_HM_equal<-list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))
freq_v_cap_HM_prop<-list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))
freq_v_cap_LM_equal<-list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))
freq_v_cap_LM_prop<-list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))

#looping over combinations, scenarios, and replicates
#saving results in 3D arrays
for(i in 1:length(combinations)) {
  #loop going through each scenario for each combination
  for(j in 1:length(scenarios)) {
    setwd(paste0(my_dir,combinations[i],scenarios[j]))
    list_files = list.files(path = paste0(my_dir,combinations[i], scenarios[j]), pattern = ".gen$")
    #loop going through every replicate for each scenario
    for(k in 1:length(list_files)) {
      #creating a temporary genind object
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #defining population boundaries by the first individual and the last individuals
      #last individual for every population as the cumulative sum of all populations (ie., last individual for pop 1 is the sum of pop 1)
      last_ind = as.numeric(cumsum(table(temp_genind@pop)))
      #first individual of every population begins at 1, then for following populations, it is the last individual (cumulative sum) + 1
      #for example, if the last individual for pop 1 is 30, the first individual for pop 2 would be 31
      first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
      #selecting the first 5 values since we have 5 populations
      first_ind = first_ind[1:5]
      
      #defining how many individuals to sample from every population
      #it varies based on which combination the loop is in (based on sample intensity)
      #the first two combinations are high intensity sampling, the last two are low intensity
      #if high intensity, sample 10% of each population, or 30 individuals equally
      if(i == 1) {
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.1)
        sample_size_equal = c(30,30,30,30,30)
      } else { #if low intensity, sample 5% of each population , or 15 individuals equally
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.05)
        sample_size_equal = c(15,15,15,15,15)
      }
      
      #defining individuals to randomly sample from each population for both strategies
      #from above, we defined the population boundaries of each population by the first and last individuals of the populations
      #ie., the 'range' of the populations is defined as (first_ind:last_ind)
      #here, we are defining which individuals to sample from, for each of the populations
      #so we are sampling a certain amount (define in sample_size... variable) from the population 'ranges'
      #also, setting a seed value for replicable results
      rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]), sample(first_ind[5]:last_ind[5], sample_size_equal[5]))
      rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]), sample(first_ind[5]:last_ind[5], sample_size_prop[5]))
      
      ##Start code for allelic capture 
      #First, calculate number of individuals per population
      n_ind <- table(temp_genind@pop)
      
      ##Then create a genpop file for temp_genind
      Spp_tot_genpop <- genind2genpop(temp_genind)
      
      ##separate by population 
      Spp_tot_genind_sep <- seppop(temp_genind)
      
      ###Start sampling code 
      ##determine alleles captured by the equal sampling 
      alleles_cap_equal <- colSums(temp_genind@tab[rows_to_samp_equal,], na.rm = T)
      
      ##determine number of alleles captured by the proportional sampling 
      alleles_cap_prop <- colSums(temp_genind@tab[rows_to_samp_prop,], na.rm = T)
      
      #Get every allele and its capture rate and frequency
      if (i==1){
        freq_v_cap_HM_prop[[j]]<-rbind(freq_v_cap_HM_prop[[j]],cbind(alleles_cap_prop,colSums(temp_genind@tab)/3000))
        freq_v_cap_HM_equal[[j]]<-rbind(freq_v_cap_HM_equal[[j]],cbind(alleles_cap_equal,colSums(temp_genind@tab)/3000))
      }
      #3000 is twice the population size- e.g. the number of gene copies (remember, its diploid)
      
      #Low mig
      if (i==2){
        freq_v_cap_LM_prop[[j]]<-rbind(freq_v_cap_LM_prop[[j]],cbind(alleles_cap_prop,colSums(temp_genind@tab)/3000))
        freq_v_cap_LM_equal[[j]]<-rbind(freq_v_cap_LM_equal[[j]],cbind(alleles_cap_equal,colSums(temp_genind@tab)/3000))
      }
      
      ##Now, determine the all the alleles captured in each category
      #First object: genpop file
      #Second object: population files (will spit out error if this is not included, irrelevant in this case)
      #Third object: # of regions, again irrelevant in this case
      #Fourth object: # of individuals per population
      allele_cat_tot <- get.allele.cat(Spp_tot_genpop, c(1:5), 2, n_ind)
      
      ##if high migration combination 
      if(i == 1) {
        
        ##calculate the total number of alleles in each frequency category over 9 allele categories
        ##array dims:
        #dim 1: replicates of each scenario (100 total)
        #dim 2: allele frequency category of all alleles - global, globally v common, etc. 
        #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
        for (a in 1:length(allele_cat_tot)) highmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
        
        ##calculate the number of alleles captured in each category for equal and proportional sampling strategies 
        ##array dims:
        #dim 1: replicates of each scenario (100 total)
        #dim 2: number of alleles captured in each category (global, globally v common, etc.) 
        #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
        for (b in 1:length(allele_cat_tot)) highmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
        for (c in 1:length(allele_cat_tot)) highmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
        
        ##calculation of the high migration percent alleles captured per category in either equal or proportional sampling strategy
        ##array dims:
        #dim 1: replicates of each scenario (100 total)
        #dim 2: % of alleles captured in each category (global, globally v common, etc.) 
        #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
        for (d in 1:length(allele_cat_tot)) highmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
        for (e in 1:length(allele_cat_tot)) highmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
        
      } #else { #if low migration
      ##calculate the total number of alleles in each frequency category over 9 allele categories for the low migration combination
      ##array dims:
      #dim 1: replicates of each scenario (100 total)
      #dim 2: allele frequency category of all alleles - global, globally v common, etc. 
      #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
      for (a in 1:length(allele_cat_tot)) lowmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
      
      ##calculate the number of alleles captured in each category for equal and proportional sampling strategies 
      ##array dims:
      #dim 1: replicates of each scenario (100 total)
      #dim 2: number of alleles captured in each category (global, globally v common, etc.) 
      #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
      for (b in 1:length(allele_cat_tot)) lowmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
      for (c in 1:length(allele_cat_tot)) lowmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
      
      ##calculation of the high migration percent alleles captured per category in either equal or proportional sampling strategy
      ##array dims:
      #dim 1: replicates of each scenario (100 total)
      #dim 2: % of alleles captured in each category (global, globally v common, etc.) 
      #dim 3: scenarios 1 - 9 (very different population sizes --> equal population sizes)
      for (d in 1:length(allele_cat_tot)) lowmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
      for (e in 1:length(allele_cat_tot)) lowmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
      
    }
  }
}

#remove row 1 which was just 0,0 
for (k in 1:length(freq_v_cap_LM_equal)) {freq_v_cap_LM_equal[[k]]<-freq_v_cap_LM_equal[[k]][-1,]; freq_v_cap_LM_prop[[k]]<-freq_v_cap_LM_prop[[k]][-1,]}
for (k in 1:length(freq_v_cap_HM_equal)) {freq_v_cap_HM_equal[[k]]<-freq_v_cap_HM_equal[[k]][-1,]; freq_v_cap_HM_prop[[k]]<-freq_v_cap_HM_prop[[k]][-1,]}
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims_bottlenecks\\bottleneck_2\\R-scripts")
save(freq_v_cap_LM_equal,freq_v_cap_LM_prop,freq_v_cap_HM_equal,freq_v_cap_HM_prop,file="Rosenberger_freq_v_cap.Rdata")

#################################################
######## High migration data frames #############
#################################################

##Loop to calculate means over all replicates for allele capture code
#store results in multiple 2D arrays
for(j in 1:length(allele_cat_tot)) {
  for (f in 1:length(allele_cat_tot)) {
    
    ##Calculate, for the high mig combo, the mean number of alleles per category (from 100 replicates)
    highmig_alleles_existing_by_cat[f,j] <-  round(mean(highmig_all_existing_by_sp_reps[,j,f]), digits = 4)
    
    ##Calculate the number of alleles in each category captured by each strategy
    highmig_equal_mean_all_cap[f,j] <- round(mean(highmig_equal_alleles[,j,f]), 4)
    highmig_prop_mean_all_cap[f,j] <- round(mean(highmig_prop_alleles[,j,f]), 4)
    
    ##Calculate the percent capture of existing alleles in these categories 
    highmig_equal_mean_all_cap_per[f,j] <- round(mean(highmig_equal_per[,j,f]), 4)*100
    highmig_prop_mean_all_cap_per[f,j] <- round(mean(highmig_prop_per[,j,f]), 4)*100 
    
  }
}

##name rows and columns 
rownames(highmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                               "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(highmig_alleles_existing_by_cat) <- list_allele_cat

setwd("C:\\Users\\kayle\\Documents\\XXX-XXX\\samp_pop_sims\\R scripts")
write.csv(highmig_alleles_existing_by_cat, "highmig_alleles_existing_by_cat.csv")

###Create data frames with percent and # of alleles captured per category

##high mig, allelic capture equal sampling strategy data frame 
highmig_all_cap_equal_df <- matrix(nrow = length(allele_cat_tot), 
                                   ncol = length(allele_cat_tot))

##loop to merge information
for(m in 1:length(highmig_equal_mean_all_cap_per[,1])){
  for(n in 1:length(highmig_equal_mean_all_cap_per[,1])){
    highmig_all_cap_equal_df[m,n] <- paste0(highmig_equal_mean_all_cap_per[m,n], "%", " ", "(", highmig_equal_mean_all_cap[m,n], ")")
  }
}  

##name rows and columns
rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                        "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_equal_df) <- list_allele_cat

##high mig, allelic capture proportional sampling strategy data frame
highmig_all_cap_prop_df <- matrix(nrow = length(allele_cat_tot), 
                                  ncol = length(allele_cat_tot))

##loop to merge information
for(o in 1:length(highmig_prop_mean_all_cap_per[,1])){
  for(p in 1:length(highmig_prop_mean_all_cap_per[,1])){
    highmig_all_cap_prop_df[o,p] <- paste0(highmig_prop_mean_all_cap_per[o,p], "%", " ", "(", highmig_prop_mean_all_cap[o,p], ")")
  }
}  

##name rows and columns
rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                        "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_prop_df) <- list_allele_cat

setwd("C:\\Users\\kayle\\Documents\\XXX-XXX\\samp_pop_sims\\R scripts")
##write out data frames
write.csv(highmig_all_cap_equal_df, "highmig_all_cap_equal_df.csv")
write.csv(highmig_all_cap_prop_df, "highmig_all_cap_prop_df.csv")

#################################################
######### Low migration data frames #############
#################################################

##Loop to calculate means over all replicates for allele capture code
#store results in multiple 2D arrays
for(j in 1:length(allele_cat_tot)) {
  for (f in 1:length(allele_cat_tot)) {
    
    ##Calculate, for the high mig combo, the mean number of alleles per category (from 100 replicates)
    lowmig_alleles_existing_by_cat[f,j] <-  round(mean(lowmig_all_existing_by_sp_reps[,j,f]), 4)
    
    ##Calculate the number of alleles in each category captured by each strategy
    lowmig_equal_mean_all_cap[f,j] <- round(mean(lowmig_equal_alleles[,j,f]), 4)
    lowmig_prop_mean_all_cap[f,j] <- round(mean(lowmig_prop_alleles[,j,f]), 4)
    
    ##Calculate the percent capture of existing alleles in these categories 
    lowmig_equal_mean_all_cap_per[f,j] <- round(mean(lowmig_equal_per[,j,f]), 4)*100
    lowmig_prop_mean_all_cap_per[f,j] <- round(mean(lowmig_prop_per[,j,f]), 4)*100 
    
  }
}

##name rows and columns 
rownames(lowmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(lowmig_alleles_existing_by_cat) <- list_allele_cat


####migration data frames to combine # and % of allelic frequencies captured 
##low mig equal sampling strategy allelic capture data frame
lowmig_all_cap_equal_df <- matrix(nrow = length(allele_cat_tot), 
                                  ncol = length(allele_cat_tot))

##loop to merge information
for(m in 1:length(lowmig_equal_mean_all_cap_per[,1])){
  for(n in 1:length(lowmig_equal_mean_all_cap_per[,1])){
    lowmig_all_cap_equal_df[m,n] <- paste0(lowmig_equal_mean_all_cap_per[m,n], "%", " ", "(", lowmig_equal_mean_all_cap[m,n], ")")
  }
}  

##name rows and columns
rownames(lowmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                       "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_equal_df) <- list_allele_cat

##low mig proportional sampling strategy allelic capture data frame
lowmig_all_cap_prop_df <- matrix(nrow = length(allele_cat_tot), 
                                 ncol = length(allele_cat_tot))

##loop to merge information
for(o in 1:length(lowmig_prop_mean_all_cap_per[,1])){
  for(p in 1:length(lowmig_prop_mean_all_cap_per[,1])){
    lowmig_all_cap_prop_df[o,p] <- paste0(lowmig_prop_mean_all_cap_per[o,p], "%", " ", "(", lowmig_prop_mean_all_cap[o,p], ")")
  }
}  

rownames(lowmig_all_cap_prop_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                      "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_prop_df) <- list_allele_cat

##write out data frames
setwd("C:\\Users\\kayle\\Documents\\XXX-XXX\\samp_pop_sims\\R scripts")
write.csv(lowmig_all_cap_equal_df, "lowmig_all_cap_equal_df.csv")
write.csv(lowmig_all_cap_prop_df, "lowmig_all_cap_prop_df.csv")

