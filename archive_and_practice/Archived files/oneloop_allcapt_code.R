#Code written by Emily Schumacher and Dr. Sean Hoban

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

###################################
###### Load Directory #############
###################################

#containing sub-folders
my_dir = "C:\\Users\\kayle\\Documents\\Morton-REU\\samp_pop_sims\\Simulations"
setwd(my_dir)
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

###function
source("G:/Shared drives/Emily_Schumacher/ten_oaks_gen/Fa_sample_funcs.R")
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

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig", "\\lowMig")

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

#####################################################
############## Run Sampling Code ####################
#####################################################

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



#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  #loop going through each scenario for each combination
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
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
      rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]), sample(first_ind[5]:last_ind[5], sample_size_equal[5]))
      rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]), sample(first_ind[5]:last_ind[5], sample_size_prop[5]))
      
      #defining how many alleles to sample based on rows_to_samp
      n_ind <- table(temp_genind@pop)
      
      ##create genepop file
      Spp_tot_genpop <- genind2genpop(temp_genind)
      
      ##separate by population 
      Spp_tot_genind_sep <- seppop(temp_genind)
      
      ##determine number of alleles captured by the equal sampling 
      alleles_cap_equal <- colSums(temp_genind@tab[rows_to_samp_equal,], na.rm = T)
      ##determine number of alleles captured by the equal sampling 
      alleles_cap_prop <- colSums(temp_genind@tab[rows_to_samp_prop,], na.rm = T)
      ##determine the total number of alleles in every category 
      allele_cat_tot <- get.allele.cat(Spp_tot_genpop, c(1:5), 2, n_ind)
      
      #if high migration 
      if(i == 1) {
        #saving number of alleles total in high migration scenario
        for (a in 1:length(allele_cat_tot)) highmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
        
        ##all the high migration allele # per category
        for (b in 1:length(allele_cat_tot)) highmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
        for (c in 1:length(allele_cat_tot)) highmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
       
         ##high migration percent alleles captured per category
        for (d in 1:length(allele_cat_tot)) highmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
        for (e in 1:length(allele_cat_tot)) highmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
      
  } #else { #if low migration
    #saving number of alleles total in low migration scenario
    for (a in 1:length(allele_cat_tot)) lowmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
    
    ##saving # of alleles captured per category   
    for (b in 1:length(allele_cat_tot)) lowmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
    for (c in 1:length(allele_cat_tot)) lowmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
    ##saving percent  of alleles per category
    for (d in 1:length(allele_cat_tot)) lowmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
    for (e in 1:length(allele_cat_tot)) lowmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
    
      }
    }
  }
#}

######################################################
### write out data files 
setwd("G:\\Shared drives\\Emily_Schumacher\\simulation_code\\allele_capture_df")

##existing alleles by species 
for(j in 1:length(highmig_alleles_existing_by_cat[1,])) {
  for (f in 1:length(allele_cat_tot)) {

      highmig_alleles_existing_by_cat[f,j] <-  round(mean(highmig_all_existing_by_sp_reps[,j,f]), digits = 4)
      highmig_equal_mean_all_cap[f,j] <- round(mean(highmig_equal_alleles[,j,f]), 4)
      highmig_prop_mean_all_cap[f,j] <- round(mean(highmig_prop_alleles[,j,f]), 4)

##calculate percent of allele capture
     highmig_equal_mean_all_cap_per[f,j] <- round(mean(highmig_equal_per[,j,f]), 4)*100
     highmig_prop_mean_all_cap_per[f,j] <- round(mean(highmig_prop_per[,j,f]), 4)*100 

  }
}

##name rows and columns 
rownames(highmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                               "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(highmig_alleles_existing_by_cat) <- list_allele_cat

write.csv(highmig_alleles_existing_by_cat, "highmig_alleles_existing_by_cat2.csv")

##create a new df = equal sampling
highmig_all_cap_equal_df <- matrix(nrow = length(highmig_equal_mean_all_cap_per[,1]), 
        ncol = length(highmig_equal_mean_all_cap_per[,1]))

  for(m in 1:length(highmig_equal_mean_all_cap_per[,1])){
    for(n in 1:length(highmig_equal_mean_all_cap_per[,1])){
       highmig_all_cap_equal_df[m,n] <- paste0(highmig_equal_mean_all_cap_per[m,n], "%", " ", "(", highmig_equal_mean_all_cap[m,n], ")")
     }
   }  

rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
   "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_equal_df) <- list_allele_cat

##make new table using proportional sampling 
highmig_all_cap_prop_df <- matrix(nrow = length(highmig_prop_mean_all_cap_per[,1]), 
       ncol = length(highmig_prop_mean_all_cap_per[,1]))


 for(o in 1:length(highmig_prop_mean_all_cap_per[,1])){
   for(p in 1:length(highmig_prop_mean_all_cap_per[,1])){
     highmig_all_cap_prop_df[o,p] <- paste0(highmig_prop_mean_all_cap_per[o,p], "%", " ", "(", highmig_prop_mean_all_cap[o,p], ")")
   }
 }  

rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
   "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_prop_df) <- list_allele_cat
  
##write out data frames
write.csv(highmig_all_cap_equal_df, "highmig_all_cap_equal_df2.csv")
write.csv(highmig_all_cap_prop_df, "highmig_all_cap_prop_df2.csv")

################################
###### Low migration tables ####
################################

##existing alleles by species 
for(j in 1:length(lowmig_alleles_existing_by_cat[1,])) {
  for (f in 1:length(allele_cat_tot)) {

    lowmig_alleles_existing_by_cat[f,j] <-  round(mean(lowmig_all_existing_by_sp_reps[,j,f]), digits = 4)
    lowmig_equal_mean_all_cap[f,j] <- round(mean(lowmig_equal_alleles[,j,f]), 4)
    lowmig_prop_mean_all_cap[f,j] <- round(mean(lowmig_prop_alleles[,j,f]), 4)

  ##calculate percent of allele capture
    lowmig_equal_mean_all_cap_per[f,j] <- round(mean(lowmig_equal_per[,j,f]), 4)*100
    lowmig_prop_mean_all_cap_per[f,j] <- round(mean(lowmig_prop_per[,j,f]), 4)*100 

  }
}
##name rows and columns 
rownames(lowmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                                 "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(lowmig_alleles_existing_by_cat) <- list_allele_cat

write.csv(lowmig_alleles_existing_by_cat, "lowmig_alleles_existing_by_cat2.csv")

##create a new df = equal sampling
lowmig_all_cap_equal_df <- matrix(nrow = length(lowmig_equal_mean_all_cap_per[,1]), 
                                  ncol = length(lowmig_equal_mean_all_cap_per[,1]))

for(m in 1:length(lowmig_equal_mean_all_cap_per[,1])){
   for(n in 1:length(lowmig_equal_mean_all_cap_per[,1])){
     lowmig_all_cap_equal_df[m,n] <- paste0(lowmig_equal_mean_all_cap_per[m,n], "%", " ", "(", lowmig_equal_mean_all_cap[m,n], ")")
    }
}  

rownames(lowmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                        "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_equal_df) <- list_allele_cat

##make new table using proportional sampling 
lowmig_all_cap_prop_df <- matrix(nrow = length(lowmig_prop_mean_all_cap_per[,1]), 
                                    ncol = length(lowmig_prop_mean_all_cap_per[,1]))


  for(o in 1:length(lowmig_prop_mean_all_cap_per[,1])){
   for(p in 1:length(lowmig_prop_mean_all_cap_per[,1])){
     lowmig_all_cap_prop_df[o,p] <- paste0(lowmig_prop_mean_all_cap_per[o,p], "%", " ", "(", lowmig_prop_mean_all_cap[o,p], ")")
   }
 }  

rownames(lowmig_all_cap_prop_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                         "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_prop_df) <- list_allele_cat

##write out data frames
write.csv(lowmig_all_cap_equal_df, "lowmig_all_cap_equal_df2.csv")
write.csv(lowmig_all_cap_prop_df, "lowmig_all_cap_prop_df2.csv")

