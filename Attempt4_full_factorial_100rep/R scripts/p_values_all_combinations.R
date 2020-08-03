#****************************************************************************************************************************************************************
#combining all p-values into one matrix
setwd("C:\\Users\\kayle\\Documents\\Morton-REU\\Attempt4_full_factorial_100rep\\R scripts")
load("p_values_highMig_highSamp.Rdata")
load("p_values_lowMig_highSamp.Rdata")
load("p_values_highMig_lowSamp.Rdata")
load("p_values_lowMig_lowSamp.Rdata")

p_value_matrix = matrix(0, nrow = 4, ncol = 9)
p_value_matrix[1,] = p_values_highMig_highSamp
p_value_matrix[2,] = p_values_lowMig_highSamp
p_value_matrix[3,] = p_values_highMig_lowSamp
p_value_matrix[4,] = p_values_lowMig_lowSamp

rownames(p_value_matrix) = c("highMig_highSamp", "lowMig_highSamp", "highMig_lowSamp", "lowMig_lowSamp")
colnames(p_value_matrix) = c("scenario 1", "scenario 2", "scenario 3", "scenario 4", "scenario 5", "scenario 6", "scenario 7", "scenario 8", "scenario 9")
round(p_value_matrix, 5)

#*****************************************************************************************************************************************************************
#loading in total alleles
load("total_alleles_highMig_highSamp.Rdata")
load("total_alleles_lowMig_highSamp.Rdata")
load("total_alleles_highMig_lowSamp.Rdata")
load("total_alleles_lowMig_lowSamp.Rdata")

mean(total_alleles_highMig_highSamp_equal)
mean(total_alleles_lowMig_highSamp_equal)

#****************************************************************************************************************************************************************
#expected heterozygosity
load("heterozygosity_highMig_highSamp.Rdata")
load("heterozygosity_lowMig_highSamp.Rdata")

heterozygosity_highMig_lowSamp_equal
heterozygosity_lowMig_lowSamp_equal
