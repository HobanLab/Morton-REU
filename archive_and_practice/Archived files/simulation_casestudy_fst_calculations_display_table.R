#Code written by Emily Schumacher and Dr. Sean Hoban in collaboration

########################################
########### Libraries ##################
########################################

library(adegenet)
library(diveRsity)
library(hierfstat)

######################################
########## Calculate PW Fst ##########
######################################
quac_dir <- 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_acerifolia'
quen_dir <- 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_engelmannii'
quog_dir <- 'C:\\Users\\kayle\\Documents\\Morton-REU\\case_study_sims\\Simulations\\q_oglethorpensis'

species_names <- c("q_acerifolia", "q_engelmannii", "q_oglethorpensis")

##final df with combined species
species_pwfst_df <- matrix(nrow = 3, ncol = 3)

##########QUAC
setwd(quac_dir)
quac_list <- list.files(quac_dir,".gen$")

quac_genind_list <- list() 

quac_hierfstat <- list()

##quac pwfst array
quac_pwfst_array <- array(dim = c(4,4,100))

##min, max, mean of replicates 
quac_mean_max_min_fst <- matrix(nrow = 3, ncol = 100)

#loop going through every replicate for each scenario
for(k in 1:length(quac_list)) {
  
  ##creating genind list for QUAC genind 
  quac_genind_list[[k]] <- read.genepop(quac_list[[k]], ncode=3)
  
  ##convert genind files to hierfstat format to run pwfst 
  quac_hierfstat[[k]] <- genind2hierfstat(quac_genind_list[[k]])
  
  ##array to store all pwfst values
  quac_pwfst_array[,,k] <- pairwise.neifst(quac_hierfstat[[k]])
  
  ##calculate statistics for QUAC - max, min, mean fst 
  quac_mean_max_min_fst[1,k] <- mean(quac_pwfst_array[,,k], na.rm = TRUE)
  quac_mean_max_min_fst[2,k] <- min(quac_pwfst_array[,,k], na.rm = TRUE)
  quac_mean_max_min_fst[3,k] <- max(quac_pwfst_array[,,k], na.rm = TRUE)
  
}
    
################QUEN
setwd(quen_dir)
quen_list <- list.files(quen_dir,".gen$")

quen_genind_list <- list() 

quen_hierfstat <- list()

##quen pwfst array
quen_pwfst_array <- array(dim = c(4,4,100))

##min, max, mean of replicates 
quen_mean_max_min_fst <- matrix(nrow = 3, ncol = 100)

#loop going through every replicate for each scenario
for(k in 1:length(quen_list)) {
  
  ##creating genind list for QUEN genind 
  quen_genind_list[[k]] <- read.genepop(quen_list[[k]], ncode=3)
  
  ##convert genind files to hierfstat format to run pwfst 
  quen_hierfstat[[k]] <- genind2hierfstat(quen_genind_list[[k]])
  
  ##calculate statistics for QUEN - max, min, mean fst 
  quen_pwfst_array[,,k] <- pairwise.neifst(quen_hierfstat[[k]])
  
  ##calculate statistics for QUEN
  quen_mean_max_min_fst[1,k] <- mean(quen_pwfst_array[,,k], na.rm = TRUE)
  quen_mean_max_min_fst[2,k] <- min(quen_pwfst_array[,,k], na.rm = TRUE)
  quen_mean_max_min_fst[3,k] <- max(quen_pwfst_array[,,k], na.rm = TRUE)
  
}

################QUOG
setwd(quog_dir)
quog_list <- list.files(quog_dir,".gen$")

quog_genind_list <- list() 

quog_hierfstat <- list()

##quog pwfst array
quog_pwfst_array <- array(dim = c(5,5,100))

##min, max, mean of replicates 
quog_mean_max_min_fst <- matrix(nrow = 3, ncol = 100)

#loop going through every replicate for each scenario
for(k in 1:length(quog_list)) {
  
  ##creating genind list for QUOG genind 
  quog_genind_list[[k]] <- read.genepop(quog_list[[k]], ncode=3)
  
  ##convert genind files to hierfstat format for pwfst
  quog_hierfstat[[k]] <- genind2hierfstat(quog_genind_list[[k]])
  
  ##QUOG pwfst array 
  quog_pwfst_array[,,k] <- pairwise.neifst(quog_hierfstat[[k]])
  
  ##calculate statistics for QUOG
  quog_mean_max_min_fst[1,k] <- mean(quog_pwfst_array[,,k], na.rm = TRUE)
  quog_mean_max_min_fst[2,k] <- min(quog_pwfst_array[,,k], na.rm = TRUE)
  quog_mean_max_min_fst[3,k] <- max(quog_pwfst_array[,,k], na.rm = TRUE)
  
}

###################################################################################

###write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
for(a in 1:length(quog_mean_max_min_fst[,1])){
  species_pwfst_df[1,a] <-  round(mean(quac_mean_max_min_fst[a,]),3)  
  species_pwfst_df[2,a] <-  round(mean(quen_mean_max_min_fst[a,]),3)
  species_pwfst_df[3,a] <- round(mean(quog_mean_max_min_fst[a,]),3)
}

##rename dfs 
rownames(species_pwfst_df) <- c("QUAC","QUEN","QUOG")
colnames(species_pwfst_df) <- c("MeanFst", "MinFst", "MaxFst")

##write out 
write.csv(species_pwfst_df,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\casestudy_pwfst_df.csv")

