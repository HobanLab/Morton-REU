#combination indicates high mig/low and high intensity samp/low
#################
#you might not necessarily need this outer loop, and it might depend on how your files are organized
# i had four files for each combination of migration rate and sample intensity 
#but this could be like 4 separater loops with just the two lower loops below (make a separate loop for high mig/high samp, low mig/low samp, etc...)
for(i in 1:length(combinations)) { 
  #####################
  #as an example -- combination 1 -> high migration and high sample intensity
  
  #inner loop going through each scenario for each combination
  for(j in 1:length(scenarios)) { 
    ##example -- scenario 1
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    #listing all .gen files to then be converted to genind objects
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    
    #convert to genind objects here and save as a list, as you already did in your code in the beginning for loops -- i called it list_genind
    
    #loop going through every genind object for each scenario (like your line 230)
    for(k in 1:length(list_genind)) {
      
      #defining population boundaries by the first individual and the last individuals 
      ######(i think you did this by doing genpop -> seppop so might not need this part)
      last_ind = as.numeric(cumsum(table(temp_genind@pop)))
      first_ind = as.numeric(c(1, cumsum(table(temp_genind@pop)) +1))
      first_ind = first_ind[1:5]
      
      #defining how many individuals to sample from every population
      ###you already did this with the rows_to_samp code, this was a more condensed way I tried doing it
      ##this relies on the outer-most loop I wrote-- for the different 'combinations' so if you don't have your files set up like that it might be different
      #high intensity -- combinations 1 and 2 were high intensity sampling
      if((i == 1) || (i == 2)) {
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.1)
        sample_size_equal = c(30,30,30,30,30)
      } else { #low intensity -- combinations 3 and 4 were low intensity sampling 
        sample_size_prop = as.numeric(table(temp_genind@pop)*0.05)
        sample_size_equal = c(15,15,15,15,15)
      }
      
      
      ###########
      #my sampling code -- again a substitute for rows_to_samp
      #defining individuals to randomly sample from each population for both strategies
      rows_to_samp_equal = c(sample(first_ind[1]:last_ind[1], sample_size_equal[1]), sample(first_ind[2]:last_ind[2], sample_size_equal[2]), sample(first_ind[3]:last_ind[3], sample_size_equal[3]), sample(first_ind[4]:last_ind[4], sample_size_equal[4]), sample(first_ind[5]:last_ind[5], sample_size_equal[5]))
      rows_to_samp_prop = c(sample(first_ind[1]:last_ind[1], sample_size_prop[1]), sample(first_ind[2]:last_ind[2], sample_size_prop[2]), sample(first_ind[3]:last_ind[3], sample_size_prop[3]), sample(first_ind[4]:last_ind[4], sample_size_prop[4]), sample(first_ind[5]:last_ind[5], sample_size_prop[5]))
      
      sample_n_alleles_equal = sum(colSums(temp_genind@tab[rows_to_samp_equal,])>0)
      sample_n_alleles_prop = sum(colSums(temp_genind@tab[rows_to_samp_prop,])>0)
      
      #saving the total alleles present in a simple variable
      total_alleles = ncol(temp_genind@tab)
      
      #####################
      #SAVING RESULTS
      #this was how I saved my results in a matrix, its a bit confusing since I have 3 nested loops
      # variable 'i' indicates the combination
      #for each combination I have separate matrices for equal and proportional strategies which I later merged, but you may not need to
      if(i == 1) {
        #saving proportion of alleles captured for both equal and proportional strategies
        ##############
        #My matrices are jxk dimensions (j = # of scenarios, so 9, and k = # of replicates, or in this edited case, # of genind objects, so 100) -- so 9X100
        #each result is saved separately for each replicate or each genind object 
        results_highMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles 
        results_highMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_highMig[j,k] = total_alleles
        
        #saving expected heterozygosity
        sum_temp_genind = summary(temp_genind)
        hexp_highMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 2) { #if low migration
        #saving proportion of alleles captured 
        results_lowMig_highSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_highSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
        #saving total alleles
        total_alleles_lowMig[j,k] = total_alleles
        
        #saving expected heterozyogisty
        sum_temp_genind = summary(temp_genind)
        hexp_lowMig[j,k] = mean(sum_temp_genind$Hexp)
        
      } else if (i == 3) {
        #saving prop alleles
        results_highMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_highMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
        
      } else {
        #saving prop alleles
        results_lowMig_lowSamp_equal[j,k] = sample_n_alleles_equal/total_alleles
        results_lowMig_lowSamp_prop[j,k] = sample_n_alleles_prop/total_alleles
      }
    }
  }
}