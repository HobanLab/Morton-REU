#---------------------------------------
#THIS IS ALL NOW IN THE SAMPLING SCRIPT, TO DELETE SOON#

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

#Place this below line 158, which is where alleles_cap_prop and alleles_cap_equal are created

###################
#High Mig
###################
#3000 is twice the population size- e.g. the number of gene copies (remember, its diploid)
if (i==1){
freq_v_cap_HM_prop[[j]]<-rbind(freq_v_cap_HM_prop[[j]],cbind(alleles_cap_prop,colSums(temp_genind@tab)/3000))
freq_v_cap_HM_equal[[j]]<-rbind(freq_v_cap_HM_equal[[j]],cbind(alleles_cap_equal,colSums(temp_genind@tab)/3000))
}

#Low mig
if (i==2){
freq_v_cap_LM_prop[[j]]<-rbind(freq_v_cap_LM_prop[[j]],cbind(alleles_cap_prop,colSums(temp_genind@tab)/3000))
freq_v_cap_LM_equal[[j]]<-rbind(freq_v_cap_LM_equal[[j]],cbind(alleles_cap_equal,colSums(temp_genind@tab)/3000))
}

#remove row 1 which was just 0,0 
for (k in 1:length(freq_v_cap_LM_equal)) {freq_v_cap_LM_equal[[k]]<-freq_v_cap_LM_equal[[k]][-1,]; freq_v_cap_LM_prop[[k]]<-freq_v_cap_LM_prop[[k]][-1,]}
for (k in 1:length(freq_v_cap_HM_equal)) {freq_v_cap_HM_equal[[k]]<-freq_v_cap_HM_equal[[k]][-1,]; freq_v_cap_HM_prop[[k]]<-freq_v_cap_HM_prop[[k]][-1,]}
save(freq_v_cap_LM_equal,freq_v_cap_LM_prop,freq_v_cap_HM_equal,freq_v_cap_HM_prop,file="Rosenberger_freq_v_cap.Rdata")

 #---------------------------------------
#ABOVE IS ALL NOW IN THE SAMPLING SCRIPT, TO DELETE SOON#

 

##################
# Analysis, make separate R file
##################

setwd("C:/Users/shoban/Documents/GitHub/Morton-REU/samp_pop_sims/")
load(file="Rosenberger_freq_v_cap.Rdata")

for (i in 1:2){
	if (i==1) { 
		this_freq_v_cap_prop<-freq_v_cap_HM_prop; this_freq_v_cap_equal<-freq_v_cap_HM_equal 
		file_append<-"_HM"
		}
	if (i==2) { 
		this_freq_v_cap_prop<-freq_v_cap_LM_prop; this_freq_v_cap_equal<-freq_v_cap_LM_equal 
		file_append<-"_LM"
		}
	
	#########################
	# ALLELE FREQS EXISTING
	#########################
	
	#This plot looks at the frequency of all alleles.  There are just more very rare alleles in scenario 1 (one large population)
	c1 <- rgb(0,102,204,max = 255, alpha = 80, names = "lt.blue")
	c2 <- rgb(255,40,40, max = 255, alpha = 80, names = "lt.pink")
	pdf(paste0("hist_allele_freqs_scen_1_v_9",file_append,".pdf"),height=5,width=8)
	hist(this_freq_v_cap_prop[[1]][,2],col=c1,breaks=500,xlim=c(0,0.05),freq=F,ylab="Density of occurrences", xlab="allele frequency bin (rarer alleles to the left)")
	hist(this_freq_v_cap_prop[[9]][,2],add=T,col=c2,breaks=500,freq=F)
	legend(0.03,30,c("Scenario 1","Scenario 9"),fill=c(c1,c2))
	dev.off()

	#Interesting side note, the other end of the spectrum Scenario 1 also has more very high frequency alleles, probably in the small pop'ns
	hist(this_freq_v_cap_prop[[1]][,2],col=c1,breaks=100,xlim=c(.3,.7),freq=F,ylab="Density of occurrences", xlab="allele frequency bin (rarer alleles to the right)",ylim=c(0,1))
	hist(this_freq_v_cap_prop[[9]][,2],add=T,col=c2,breaks=100,freq=F)

	#########################
	# ALLELE FREQS CAPTURED OR NOT
	#########################
	
	#This plot shows that there are more higher frequency alleles NOT captured (0 counts) for scenario 1, when equal sampling v. proportional; the plus/ minus is a neat way to show 'jitter'
	plot(this_freq_v_cap_prop[[1]][,1]-.15,this_freq_v_cap_prop[[1]][,2],ylim=c(0,.1),xlim=c(0,10))
	abline(lm(this_freq_v_cap_prop[[1]][,2]~this_freq_v_cap_prop[[1]][,1]))
	points(this_freq_v_cap_equal[[1]][,1]-.05,this_freq_v_cap_equal[[1]][,2],col="red")
	abline(lm(this_freq_v_cap_equal[[1]][,2]~this_freq_v_cap_equal[[1]][,1]),col="red")
	points(this_freq_v_cap_prop[[9]][,1]+.05,this_freq_v_cap_prop[[9]][,2],col="darkgrey") 
	abline(lm(this_freq_v_cap_prop[[9]][,2]~this_freq_v_cap_prop[[9]][,1]),col="darkgrey")
	points(this_freq_v_cap_equal[[9]][,1]+.15,this_freq_v_cap_equal[[9]][,2],col="purple") 
	abline(lm(this_freq_v_cap_equal[[9]][,2]~this_freq_v_cap_equal[[9]][,1]),col="purple")

	#This shows that in more detail how more higher frequency alleles are NOT captured
	pdf(file=paste0("boxplot_freqs_missed_alleles",file_append,".pdf"),height=8,width=5.5)
	boxplot(this_freq_v_cap_equal[[1]][(this_freq_v_cap_equal[[1]][,1])==0,2],this_freq_v_cap_equal[[9]][(this_freq_v_cap_equal[[9]][,1])==0,2],
	this_freq_v_cap_prop[[1]][(this_freq_v_cap_prop[[1]][,1])==0,2],this_freq_v_cap_prop[[9]][(this_freq_v_cap_prop[[9]][,1])==0,2],
	col=c(rep("lightblue",2),rep("darkblue",2)),names=c("S1, equal","S9, equal","S1, prop","S9, prop"),ylab="Proportion of those alleles NOT sampled")
	dev.off()
	t.test(this_freq_v_cap_equal[[1]][(this_freq_v_cap_equal[[1]][,1])==0,2],this_freq_v_cap_equal[[9]][(this_freq_v_cap_equal[[9]][,1])==0,2])
	 # mean of x   mean of y 
	#0.003172394 0.002118881 		#Almost 50% higher
	
	#  mean of x   mean of y 
	#0.009902135 0.004192742 
	hist(this_freq_v_cap_equal[[9]][(this_freq_v_cap_equal[[9]][,1])==0,],breaks=50,col=c2,freq=F)
	hist(this_freq_v_cap_equal[[1]][(this_freq_v_cap_equal[[1]][,1])==0,],breaks=50,col=c1,add=T,freq=F)
}