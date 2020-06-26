library(adegenet)
library(diveRsity)

setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop")

#import functions
import_arp2gen_files = function(mypath, mypattern) {
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

import_gen2genind_objects = function(mypath, mypattern) {
  temp_list_3 = list.files(mypath, mypattern)
  temp_list_4 = list(length = length(temp_list_3))
  for(j in 1:length(temp_list_3)){temp_list_4[[j]]=read.genepop(temp_list_3[j], ncode=3)}
  temp_list_4
}


##importing files
import_arp2gen_files("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop", ".arp$")

genind_list = list()
genind_list = import_gen2genind_objects("C:\\Users\\kayle\\Documents\\Morton-REU\\practice\\Simulations\\simcoal2\\example_file_2pop", ".gen$")


##analyzing/storing results
##This loop will store results in a 3D array
##calculating 2 statistics for two populations
results = array(0, dim=c(2,2,5))##note syntax
#seppop to separate each fiel into two populations
for(i in 1:length(genind_list)) {
  obj[i] = list(seppop(genind_list[[i]]))
}
for(i in 1:length(genind_list)) {
  results[1,1,i] = 1
  results[1,2,i] = 2
  results[2,1,i] = 3
  results[2,2,i] = 4
}

<<<<<<< HEAD
##plotting
for(i in 1:length(genind_list)) {
  x[i] = list(summary(genind_list[[i]]))
}

#heterozygosity expected vs. observed
plot(x[[1]]$Hexp, x[[1]]$Hobs, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)
=======
#summary
sum_rep_1 = summary(rep_1)
sum_rep_1$pop.n.all

#plotting statistics for 
plot(sum_rep_1$n.by.pop, sum_rep_1$pop.n.all, xlab="Sample size", 
     ylab="Number of alleles", main="Alleles numbers and sample sizes")
text(sum_rep_1$n.by.pop, sum_rep_1$pop.n.all, lab=names(sum_rep_1$n.by.pop))
>>>>>>> 57e945ea91795e850df2e9df66bfec28fa287316

plot(x[[2]]$Hexp, x[[2]]$Hobs, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)

plot(x[[3]]$Hexp, x[[3]]$Hobs, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)

plot(x[[4]]$Hexp, x[[4]]$Hobs, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)

plot(x[[5]]$Hexp, x[[5]]$Hobs, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
abline(0,1,lty=2)
