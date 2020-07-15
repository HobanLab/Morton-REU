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
    Figures: contains image figures (graphs in R)
Attempt3_full_factorial_lowRep: Thrid attempt for project; importing entire simulated population into R with varying parameter files (low replicates)
    R scripts: contains data analysis in R
    Simulations: contain parameter files and other simulation files
        highMig_highSamp: contains scenarios with high migration rates and high sampling intensity
        highMig_lowSamp: contains scenarios with high migration rates and low sampling intensity
        lowMig_highSamp: contains scenarios with low migration rates and high sampling intensity
        lowMig_lowSamp: contains scenarios with low migration rates and low sampling intensity
    Figures: contain image files (graphs in R)
Attempt4_full_factorial_highReps
    R scripts: contain data analysis in R (imports simulation files from outside file location)
    Figures: contains image files (graphs made in R)
    note: simulation files for this attempt are stored in an outside file location due to the larger number of replicates

Simulation files:
.par .arp .gen .simparam
Edited in text editor Notepad++
Run through the software Simcoal2
software takes .par files and converts them to .arp files (genetic data)

Rscripts:
.R
import .arp files into R for conversion to .gen files through adegenet package
Convert .gen files to genind objects through adegenet package
Analyze data through functions associated with the adegenet package
Create figures for data visualization with package ggplot2

