# Morton REU and REEF project
Simulation files and R scripts associated with summer REU project and extension fellowship at the Morton Arboretum. Code files were written by Kaylee Rosenberger, Emily Schumacher, and Dr. Sean Hoban in collaboration. 
##### Background
The overall aim of this project was to contribute to practical seed sampling guidelines for creating and maintaining genetically diverse collections for botanic garden and arboreta. Informing these sampling guidelines is one way to ensure a genetically representative sample is obtained from wild populations. Prior work has found that it is important to consider species' traits like dispersal, mode of reproduction, population history, and more, when sampling from wild populations. **For this project, we were specifically interested in creating guidelines that could apply to species with unequal population sizes.** Previous simulation studies assumed populations of equal sizes for simplicity of the model. Thus, we want to determine how to effectively capture genetic diversity when population sizes vary for a rare species. **The two 'strategies' we tested were equal (an equal number from each population regardless of size) and proportional (sampling proportionally to the population's size).**

#### Summary
There were two main aspects of this project. First, we simulated a hypothetical rare species using Simcoal 2. Then, the two sampling strategies were tested on simulated populations of varying sizes. We then created scripts with R that represent sampling from the simulated populations. The two sampling strategies were assessed based on the amount of genetic diversity captured (in proportion of alleles captured). The equal strategy takes a constant number of individuals from every population regardless of their sizes. The proportional strategy allocates sampling effort--it samples more from larger populations and less from smaller populations. Results were analyzed and visualized in R. 

#### Other parameters tested
We also tested different migration rates and sampling intensities for this portion of the project to determine how results would vary based on these values.

#### Case studies
Also, we created case studies to determine whether the same results would be achieved when using realistic parameter values that represent real species, in comparison to our first sets of simulations, which represented a hypothetical rare species. Our case studies were based on 3 species of Oak--Quercus acerifolia, Quercus engelmannii, and Quercus oglethorpensis. Similar to the framework described above, we created realistic parameter values for these species, simulated populations using Simcoal 2, and simulated sampling using R code to randomly select individuals. Again, we tested high and low sampling intensities to compare the diversity captured through each. 
    
Files include R scripts used for data collection, analysis, and producing plots, and text files containing paramter values used for simulation. 
As this started as an REU project, multiple attempts or trials were used to gradually build the project, increase realism and revise code.  These initial attempts are stored documented in separate folders, but are not needed to recreate the final results for the project. Each attempt builds on the previous in improving the code efficiency or adding new parameter values.

#### Simulating bottlenecks 
To increase the complexity of our simulations, we implemented bottlenecks. We simulated two types of bottlenecks for our first sets of simulations. The first type of bottleneck (bottleneck_1) simulates a population constriction of 10x for each population 5 generations past. The second bottleneck (bottleneck_2) simulates a population constriction of varying factors depending on the population size. All current populations (which vary by some degree in size) were historically all the same size (3000 individuals), and constriction happened at varying rates such that we get current populations of different sizes. 

#### File types
**Parameter files:**
    .par .txt
    Edited in text editor Notepad++
    These are input to the software Simcoal/Simcoal2 to create genetic datasets The .par signifies parameter files.  They contain information to create the genetic datasets via a coalescent simulation, including population sizes and migration rates.  The .txt files or 'batch files' contain a list of .par files to be run sim through the simulation software with one command

**Simulation output files:**
    .par .arp .gen .simparam
    Created through the software Simcoal2 after a parameter file is successfully imported into Simcoal2 and the simulation is run.  The .arp files (genetic data) are the initial dataset in Arlequin format; the .gen files are the datasets after conversion to genepop format.  The .simparam is just a mirror file of simulation parameters run.

**Rscripts:**
    .R 
    For this project, R scripts were used to import .arp files into R for conversion to .gen files through adegenet package, convert .gen files to genind objects through adegenet package, analyze data through functions associated with the adegenet package, and create figures for data visualization with package ggplot2 (see flowchart).  There are four R scripts in the simulation folder.  They are:

**samp_pop_sims and samp_pop_sims_bottlenecks R-script folder contents** (numbers indicate order in which to run files)
1_sampling_calcs.R: converts from Arlequin format to genepop format, loops to import genepop files to genind objects in R and then calculate and save basic genetic statistics (e.g. heterozygosity).
2_data_prep.R: Takes the sampling results and puts them in tidy format for graphing and analysis
3_plots_stats.R: Plots results, and performs Wilcoxon statistical tests, and p value adjustment for multiple comparisons.
4_all_freq_categories.R: Calculates alleles frequncy present in high and low migration scenarios, as well as the total alleles present 
5_extra_plots_all_freq_details.R: Various plots to visualize the frequency of alleles to determine how many rare alleles exist in scenario 1 vs. scenario 9. Also plots to determine what proportion of these rare alleles are 'missed' by the sampling strategy 
Fa_sample_funcs_source.R: Function written by Dr. Sean Hoban to calculate allele category frequencies. 

**Note:** For the case studies, the scripts listed above are combined into a single script for each species (files listed in the directory by species). The same general logic is performed in each script. These scripts can be performed in any order; however, the combined_graphics_analysis.R script must be performed last, after all species have been processed. The scripts for both case_study_sims and case_study_sims_bottlenecks are:
q_acerifolia.R: Contains file conversion functions, main sampling loop, genetic statistics, and single plot for Q. acerifolia
q_engelmannii.R: Contains file conversion functions, main sampling loop, genetic statistics, and single plot for Q. engelmannii
q_oglethorpensis.R: Contains file conversion functions, main sampling loop, genetic statistics, and single plot for Q. oglethorpensis
combined_graphics_analyses.R: This script plots all three case study species on a single plot


###### Directory contents:

    archive_and_practice: contains practice simulations as well as third, fourth, and fifth attempts for the project; importing entire simulated population into R with varying parameter files (varying rates of migration); various number of simulation replicates (note: These files are just to track progress. They do not contain complete analyses)

    samp_pop_sims: Updated code from the last iteration to analyse all combinations of migration and sampling intensities within one nested for loop for increased efficiency and ease. Also removed any hard-coded instances of rows_to_samp variable compared to Attempt5. (Sixth and final attempt for the project)
        R-scripts: R scripts to run sampling code on simulation files, convert files, save genetic statistics, run various analyses and plot results
        Simulations: contain parameter files and other simulation files; parameter files include high and low migration; simulation outputs are stored in folders.
            highMig:
            lowMig:
        Figures: contain images produced in R scripts for boxplots and other various figures

    samp_pop_sims_bottlenecks: Contains simulation files, R scripts, and figures for the 2 sets of bottleneck simulations. Bottlenecks were implemented to increase the complexity of our simulations. 
        bottleneck_1: This simulation constricts populations 10x from their historical size to their current size
        bottleneck_2: This simulation constricts historical populations of equal sizes at unequal factors to achieve current populations of unequal sizes. 

    case_study_sims: Contains R scripts and simulation files relating to three species of Oak that were used as case studies for the project. Case studies were used in comparison to a hypothetical species (in previous attempts) in order to determine if, when applied to mroe realistic parameter values, the same results would be achieved. 
        Figures: contain images produced in R scripts for boxplots
        R scripts: three separate R scripts for each case study species. Parameter values used were specified for each species separately in order to best represent them in simulation. Migration was constant in comparison to previous attempts, and different sampling intensities were used. Calculates Fst for case study species as well
            q_acerifolia
            q_engelmannii
            q_oglethorpensis
            combined_graphics_analyses.R: this script was used to merge results from all case studies into one large dataframe to be plotted together for comparison, Fst calculations for all species are also in this script
        Simulations: Contains parameter files used to simulate each of the case study species, along with simulation result files.

    case_study_sims_bottlenecks: Contains two folders for each of the bottleneck simulations. Has similar contents as described above, for case_study_sims
        bottleneck_1: This simulation constricts populations 10x from their historical size to their current size
        bottleneck_2: This simulation constricts historical populations of equal sizes at unequal factors to achieve current populations of unequal sizes. 

