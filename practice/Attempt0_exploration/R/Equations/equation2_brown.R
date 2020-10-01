p <-seq(0.05,0.001,by=-.001) #defining frequencies by creating a sequence of numbers from 0.05 to 0.001 by -0.001
s = rep(0, length(p)) #declaring S to be the same length as p, filled with values of 0

#for loop, iterates through the length of p calculating S for each value of p
for(i in 1:length(p)) {
  S[i] = 3/(log(1-p[i])) #re-declaring S with the values calcualted by the equation
}
plot(S~p)
#sample size required to capture at least one copy of a rare allele increases exponentially as the rareness of the allele increases


#question:
#I coded the following to try to declare values of p, but it didn't work. The error message stated that S and p weren't the same lengths:
#p = c(0.05, 0.025, 0.01215, 0.00625)
#S = rep(0, length(p))