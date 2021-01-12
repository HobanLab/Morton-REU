# Morton REU & REEF Project 2020
Simulation files and R scripts associated with summer REU project and extension fellowship at the Morton Arboretum. 
    The overall aim of this project was to contribute to practical seed sampling guidelines for creating and maintaining genetically diverse collections for botanic garden and arboreta. Informing these sampling guidelines is one way to ensure a genetically representative sample is obtained from wild populations. Prior work has found that it is important to consider species' traits like dispersal, mode of reproduction, population history, and more, when sampling from wild populations. For this project, we were specifically interested in creating guidelines that could apply to species with unequal population sizes. Previous simulation studies assumed populations of equal sizes for simplicity of the model. Thus, we want to determine how to effectively capture genetic diversity when population sizes vary for a rare species. The two 'strategies' we tested were equal (an equal number from each population regardless of size) and proportional (sampling proportionally to the population's size).
    There were two main aspects of this project. First, we simulated a hypothetical rare species using Simcoal 2. Then, the two sampling strategies were tested on simulated populations of varying sizes. We then created scripts with R that represent sampling from the simulated populations. The two sampling strategies were assessed based on the amount of genetic diversity captured (in proportion of alleles captured). The equal strategy takes a constant number of individuals from every population regardless of their sizes. The proportional strategy allocates sampling effort--it samples more from larger populations and less from smaller populations. Results were analyzed and visualized in R. We also tested different migration rates and sampling intensities for this portion of the project to determine how results would vary based on these values.
    Then, we created case studies to determine whether the same results would be achieved when using realistic parameter values that represent real species, in comparison to our first sets of simulations, which represented a hypothetical rare species. Our case studies were based on 3 species of Oak--Quercus acerifolia, Quercus engelmannii, and Quercus oglethorpensis. Similar to the framework described above, we created realistic parameter values for these species, simulated populations using Simcoal 2, and simulated sampling using R code to randomly select individuals. Again, we tested high and low sampling intensities to compare the diversity captured through each. 
    
    Files include R scripts used for data collection, analysis, and producing plots, and text files containing paramter values used for simulation. 
    Multiple attemps are documented in separate folders. Each attempt builds on the previous in improving the code efficiency or adding new parameter values.

Parameter files:
    .par .txt
    Edited in text editor Notepad++
    Input to the software Simcoal/Simcoal2 to create genetic datasets
    .par represent parameter files containing information to create the genetic datasets
    .txt files or 'batch files' contain a list of .par files to be sim through the simulation software with one command

Simulation files:
    .par .arp .gen .simparam
    Created through the software Simcoal2 after a parameter file is successfully imported into Simcoal2 and converted
    software takes .par files and converts them to .arp files (genetic data)

Rscripts:
    .R
    For this project, R scripts were used to import .arp files into R for conversion to .gen files through adegenet package, convert .gen files to genind objects through adegenet package, analyze data through functions associated with the adegenet package, and create figures for data visualization with package ggplot2

Flow chart describing order in which to run files for hypothetical species (samp_pop_sims)
![Alt text](samp_pop_sims/Figures/read_me_flowchart.png?raw=true "Files to run")
Flow chart describing order in which to run files for case study species (case_study_sims)
![Alt text](case_study_sims/Figures/read_me_flowchart_case_studies.png?raw=true "Files to run")

Directory contents:

    practice_sims: Contains practice simulations and R scripts
        Attempt0_exploration: Exploring functions and objects related to this project
            R: Following instructions on the genind tutorial ( )
            Simulations: Contains practice parameter files for the software Simcoal and Simcoal 2
        Attempt1_practice_sims: (archived to practice folder) Contains first iteration of writing my own parameter files, running simulations, and analyzing the code
            Simulations: contains parameter files used for simulation, takes a sample of the entire population
            R scripts: contains analysis files. Wrote functions to import files from a directory and convert them. Wrote a for loop to convert imported files to genind objects
        Attempt2_practice_sims_revised: Contains similar parameter files and analyses as Attempt1, but code and parameter values have been revised
            Simulations: contains parameter files used for simulation, takes a sample of the entire population, revised parameter values
            R scripts: contains analysis files. Wrote functions to import files from a directory and convert them. Wrote a for loop to convert imported files to genind objects

    archive_previous_versions: contains third, fourth, and fifth attempts for the project; importing entire simulated population into R with varying parameter files (varying rates of migration); various number of simulation replicates
        Attempt3_mig_intensity_10_reps: testing code and parameters using low level of replicates for faster loading
            Figures: contain figures generated from R
            R scripts: contains R scripts for sampling code and Rdata files for faster loading
            Simulations: contains all parameter files used for simulation, along with the results of simulation
        Attempt4_mig_intensity_100_reps: tesing code and parameters with higher replicates in order to view more accurate/precise results
            Figures: contain figures generated from R
            R scripts: contains R scripts for sampling code and Rdata files for faster loading
            Simulations: contains all parameter files used for simulation, along with the results of simulation
        Attempt5_REU_final: Finalized version of code for Summer REU project. Includes 100 simulation replicates, high and low migration rates, and high and low sampling intensity. Slightly improved code efficiency. 
            Figures: contain figures generated from R
            R scripts: contains R scripts for sampling code and Rdata files for faster loading
            Simulations: contains all parameter files used for simulation, along with the results of simulation


    samp_pop_sims: Updated code from the last iteration to analyse all combinations of migration and sampling intensities within one nested for loop for increased efficiency and ease. Also removed any hard-coded instances of rows_to_samp variable compared to Attempt5. (Sixth and final attempt for the project)
        R scripts: Improved code to analyse all combinations of migration and sampling intensities within one loop. Removed most instances of hard-coding to improve flexibility
            analysis_conversions_calculations_revised.R: revised vresion fixes errors. Data was being analyzed incorrectly.
            analysis_data_prep.R: Script to aggregate dataframes from conversion file into one main dataframe for generating results. This is to get results for equal and proportional sampling on one plot
            analysis_graphics_results.R: Script creates boxplots for each combination of migration rates and sample intensity. 
        Simulations: contain parameter files and other simulation files; parameter files include high and low migration; simulation outputs are stored in folders.
            highMig: simulation files resulting from high migration parameters
            lowMig:  simulation files resulting from low migration parameters
        Figures: contain images produced in R scripts for boxplots

    cast_study_sims: Contains R scripts and simulation files relating to three species of Oak that were used as case studies for the project. Case studies were used in comparison to a hypothetical species (in previous attempts) in order to determine if, when applied to mroe realistic parameter values, the same results would be achieved. 
        Figures: contain images produced in R scripts for boxplots
        R scripts: three separate R scripts for each case study species. Parameter values used were specified for each species separately in order to best represent them in simulation. Migration was constant in comparison to previous attempts, and different sampling intensities were used. 
            q_acerifolia
            q_engelmannii
            q_oglethorpensis
        Simulations: Contains parameter files used to simulate each of the case study species, along with simulation result files.
