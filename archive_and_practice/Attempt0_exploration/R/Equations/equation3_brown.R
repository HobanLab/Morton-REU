#defining variables with similar values to the table in Brown paper
p = c(0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.05)
#S = c(0, 0, 0, 0, 0, 0)
S = matrix(1:36, byrow = TRUE, nrow = 6)
r = c(1, 2, 4, 8, 10, 15)

#using equation 3 to calculate S
for(i in 1:length(p)) {
   for(j in 1:length(r)) {
    S[i,j] = (r[j] + (1.645*sqrt(r[j])) + 0.5)/(p[i])
  }
}

rownames(S) <- r
colnames(S) <- p

#change S to a data frame
S <- as.data.frame(S)

#print S with rownames and colnames defined by r and p
S

#neater output
round(S)
