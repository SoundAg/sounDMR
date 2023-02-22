# DMR Analysis Functionalization

## Authors
Tom Cairns  
Jack Colicchio  
Keerthana N Prabhu


## Purpose
This project intends to create useful, robust functions for DMR analysis. It takes
in a megaframe, cleans with the help of the programmer, then runs the
random effects model for a given formula, and finally identifies the changepoints
in the Z score of interest from the random effects model. 


## Files
*DMR_analysis_workflow*

This file contains the actual workflow for the optimized DMR analysis pipeline. 
It can be used either as a reference for analysis or as the script for analysis
depending on the experiment.


## Package
The sounDMR2 package is contained within the sounDMR2 directory. The current
package can be downloaded from github in the .tar.gz format.


## Guide
This README will not contain a comprehensive guide to how to run the workflow,
but rather an overview. The full guide for our analysis can be found on google
drive.

#### Install and Load Package
The first step is to install the package from this repository. Download the
.tar.gz file and then load into into RStudio by going to Tools -> Install Packages.
Once the package is installed it can be loaded into the working analysis by
using `library(sounDMR2)`.

#### Create Megaframe
The first step is to generate a megaframe that contains methyl bed information
for each individual combined into a single data frame.

#### Create Zoomframe
DMR analysis requires some additional columns apart from just the methyl bed information.
This function creates,
Zeroth_pos : Position adjusted to have position 0 as ATG for the target;
Gene : This can inclkude either Gene Id or Gene names for tracking during DMR analysis;
Zoom_co : This keeps track of whether the given position is within the gene or is an adaptive region.

#### Clean and Rearrange Data
This step cleans the data read into the environment and creates "long" formats
of the count of methylated and unmethylated reads as well as the percent
methylation for each individual.

#### Creating Methylation Summary Data
This section is involved in manipulating the methylation data to summarize
information based on specific variables such as the Percent methylation of each
individual for each cytosine.

#### Creating Differential Methylation Output File
This section creates a data frame where all the output data will be stored. It
creates three columns for each individual: the Z score, the change in methylation
between the treatment and the control, and the read depth for an individual.

#### Group DMR Analysis
This section is where the model is run for our DMR analysis. Here we use a
binomial mixed effects model where the user has the option to change the fixed
and random effects passed into the model depending on the goals of the experiment.

#### Changepoint Analysis
The output from the previous step gets passed into our changepoint analysis
section. Here the user picks a column upon which to run the analysis. We
recommend using one of the Z_score columns created from the fixed effects provided
to the model.


