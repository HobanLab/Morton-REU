# Morton-REU
Simulation files and R scripts associated with summer REU project and extension fellowship at the Morton Arboretum
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
    Created through the software Simcoal2
    software takes .par files and converts them to .arp files (genetic data)

Rscripts:
    .R
    import .arp files into R for conversion to .gen files through adegenet package
    Convert .gen files to genind objects through adegenet package
    Analyze data through functions associated with the adegenet package
    Create figures for data visualization with package ggplot2

![Alt text](casestudies/Figures/read_me_flowchart.png?raw=true "Files to run")

Directory contents:

    Practice: Contains practice simulations and R scripts
        Attempt0_exploration: Exploring functions and objects related to this project
            R: Following instructions on the genind tutorial ( )
            Simulations: Contains practice parameter files for the software Simcoal and Simcoal 2
        Attempt1_practice_sims: (archived to practice folder) Contains first iteration of writing my own parameter files, running simulations, and analyzing the code
            Simulations: contains parameter files used for simulation, takes a sample of the entire population
            R scripts: contains analysis files. Wrote functions to import files from a directory and convert them. Wrote a for loop to convert imported files to genind objects
        Attempt2_practice_sims_revised: Contains similar parameter files and analyses as Attempt1, but code and parameter values have been revised
            Simulations: contains parameter files used for simulation, takes a sample of the entire population, revised parameter values
            R scripts: contains analysis files. Wrote functions to import files from a directory and convert them. Wrote a for loop to convert imported files to genind objects

    Attempt3_mig_intensity_10_reps: Third attempt for project; importing entire simulated population into R with varying parameter files; using 10 simulated replicates
        R scripts: contains data analysis in R, script was updated to include sampling code. Code imports simualted files and converts them to genind objects. For loop iterates through each simulation replicate to simulate sampling from each population
            highMig_highSamp_analysis_combined.R: analysis for high migration parameter files, samples 10% from each population 
            highMig_lowSamp_analysis_combined.R: analysis for low migration parameter files, samples 10% from each population
            lowMig_highSamp_analysis_combined.R:  analysis for high migration parameter files, samples 5% from each population 
            lowMig_lowSamp_analysis_combined.R: analysis for low migration parameter files, samples 5% from each population
        Simulations: contain parameter files and other simulation files; parameter files were updated to include high and low migration;  simulation outputs are stored in folders
            highMig_highSamp: simulation files for high migration rates and high sampling intensity
            highMig_lowSamp: simulation files for high migration rates and low sampling intensity
            lowMig_highSamp: simulation files for low migration rates and high sampling intensity
            lowMig_lowSamp: simulation files for low migration rates and low sampling intensity
        Figures: contain plots resulted from analyses in R scripts

    Attempt4_mig_intensity_100_reps: increased replicates from 10 to 100 from last iteration. Improved analysis files.
        R scripts: contain data analysis in R scripts. Folder also contains .Rdata files, to save time with the larger amount of replicates, processed data may be loaded in without re-running code.  updated to create boxplots and compute p-vals
            highMig_highSamp_analysis_100reps.R: analysis for high migration parameter files, samples 10% from each population
            lowMig_highSamp_analysis_100reps.R: analysis for low migration parameter files, samples 10% from each population
            highMig_lowSamp_analysis_100reps.R: analysis for high migration parameter files, samples 5% from each population
            lowMig_lowSamp_analysis_100reps.R: analysis for low migration parameter files, samples 5% from each population
        Simulations: contain parameter files and other simulation files; parameter files were updated to include high and low migration; simulation outputs are stored in folders
            highMig_highSamp: simulation files for high migration rates and high sampling intensity
            highMig_lowSamp: simulation files for high migration rates and low sampling intensity
            lowMig_highSamp: simulation files for low migration rates and high sampling intensity
            lowMig_lowSamp: simulation files for low migration rates and low sampling intensity
        Figures: contains plots resulted from analyses in R scripts, contains figures used for describing populations.

    Attempt5_REU_final: 100 simulation replicates. Improved code to process all simulation files with one main script. Scripts are split into different files used for converting/calculating, preparing data and generating graphics.
        R scripts: Contains data analysis in R scripts. Folder also contains .Rdata files
            analysis_conversions_calculations.R: R script to process both high and low migration scenarios, as well as sample with high and low intensity within one file (for simplicity). Saves results to be imported into the data_prep script.
            analysis_data_prep.R: Script to aggregate dataframes from conversion file into one main dataframe for generating results. This is to get results for equal and proportional sampling on one plot
            analysis_graphics_results.R: Script creates boxplots for each combination of migration rates and sample intensity. 
        Simulations: contain parameter files and other simulation files; parameter files include high and low migration; simulation outputs are stored in folders. Updated from the last iteration by running the simulation only twice--one for high migration and one for low (not necessary to run four times, since sampling occurs within the code)
            highMig: simulation files resulting from high migration parameters
            lowMig: simulation files resulting from low migration parameters
        Figures: contain images produced in R scripts for boxplots

    Attempt6_complete_loop: Updated code from the last iteration to analyse all combinations of migration and sampling intensities within one nested for loop. Also removed any hard-coded instances of rows_to_samp variable compared to Attempt5. 
        R scripts: Improved code to analyse all combinations of migration and sampling intensities within one loop. Removed most instances of hard-coding to improve flexibility
            analysis_conversions_calculations_revised.R: revised vresion fixes errors. Data was being analyzed incorrectly.
            analysis_data_prep.R: Script to aggregate dataframes from conversion file into one main dataframe for generating results. This is to get results for equal and proportional sampling on one plot
            analysis_graphics_results.R: Script creates boxplots for each combination of migration rates and sample intensity. 
        Simulations: contain parameter files and other simulation files; parameter files include high and low migration; simulation outputs are stored in folders.
            highMig: simulation files resulting from high migration parameters
            lowMig:  simulation files resulting from low migration parameters
        Figures: contain images produced in R scripts for boxplots
