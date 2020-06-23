library(diveRsity)
library(adegenet)

##try to load in arlequin file and then convert it
setwd('C:\\Users\\kayle\\Downloads')
##convert arlequin file to a genepop format
butternut_genepop <- arp2gen("butternut_arl.arp")
##now read in gene pop file into genind for adegenet 
butternut_genind <- read.genepop("butternut_arl.gen", ncode = 3)



#examining the object
butternut_genind
#accessing slots of the object
#@tab gives a matrix of the allele counts
head(butternut_genind@tab)
#@ploidy gives ploidy for each individual
butternut_genind@ploidy

#using accessors
nInd(butternut_genind) #1761 individuals
nLoc(butternut_genind) #11 microsatellite loci
nPop(butternut_genind) #44 populations
nAll(butternut_genind) #alleles at each locus
tab(butternut_genind) #same as butternut_genind@tab ?

#estimating inbreeding
#inbreeding estimate for each individual
temp = inbreeding(butternut_genind)
##calculating the mean inbreeding
Fbar = sapply(temp, mean)
#histogram of Fbar
hist(Fbar, main="Average inbreeding")

##subsetting individuals with inbreeding estimates over 0.6
which(Fbar > 0.6)
F_ind = inbreeding(butternut_genind, res.type="function")[which(Fbar>0.6)]
F_ind
plot(F_ind$"0118", main="Inbreeding",
     xlab="Inbreeding (F_ind)", ylab="Probability density")






##convert to genpop
butternut_genpop <- genind2genpop(butternut_genind)
butternut_genpop

#matrix of values
butternut_genpop@tab

#subset first 5 populations
subset_butternut = butternut_genpop[1:5]
#checking
popNames(butternut_genpop)
popNames(subset_butternut)

#loci names
locNames(butternut_genpop)
#subsetting based on loci
subset_butternut_loc = butternut_genpop[loc = c("locus1", "locus7", "locus10", "locus11")]
#checking
locNames(subset_butternut_loc)

#trying seploc
seploc_butternut = seploc(butternut_genpop)
seploc_butternut
names(seploc_butternut)
#new object that selects just one locus
loc7_butternut = seploc_butternut$locus7
loc7_butternut

#seppop
seppop_butternut = seppop(butternut_genind)
#separates genind object by population 1
seppop_butternut$pop1

#utility of lapply
seppop_butternut = lapply(seppop_butternut, seploc)
names(seppop_butternut)
