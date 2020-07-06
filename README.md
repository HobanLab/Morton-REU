# Morton-REU
Simulation files and Rscripts associated with summer REU project

Directory:
Practice folder: practice simulations and R scripts
Attempt1_samppop: First attempt for project goals; importing simulated samples from simcoal2 into R
    R: contains R scripts
    Simulations: contains simulation files 
Attempt2_fullpop: Second attempt for project; importing the entire simulated population into R
    R scripts: contains data analysis in R
    Simulations: contains parameter files and other simulation files

Simulation files:
.par .arp .gen .simparam
Edited in text editor 
Run through the software Simcoal2
software takes .par files and converts them to .arp files (genetic data)

Rscripts:
.R
import .arp files into R for conversion to .gen files through adegenet package
Convert .gen files to genind objects through adegenet package
Analyze data through functions associated with the adegenet package

