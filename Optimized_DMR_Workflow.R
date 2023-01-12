#!/usr/bin/env

#-------------------------------------------------------------------------------

# Name: Optimized_DMR
# Author: Tom Cairns
# Purpose: This script is the workflow for the DMR analysis.

#-------------------------------------------------------------------------------

# library(reshape2)
# library(lme4)
# library(changepoint)
# library(tidyverse)
# library(Formula)
library(sounDMR2)

#-------------------------------------Main--------------------------------------

#--------------
# Read In Data
#--------------
#Load in megaframe
ZoomFrame <- read.csv(file.choose(), head=TRUE)
#Load in ID file
Exp_ID <- read.csv(file.choose(), head=TRUE)

#------------
# Clean Data
#------------

# The colnames of interest. This will be carried through the rest of the analysis
colnames_of_interest <- c('Chromosome', 'Gene', 'Position', 'Strand', 'CX', 
                          'Zeroth_pos', 'Plant')

cleaned_data <- clean_data(ZoomFrame, Exp_ID)

# Extract the cleaned data
ZoomFrame_filtered <- cleaned_data$ZoomFrame_filtered
LongPercent <- cleaned_data$LongPercent
LongMeth <- cleaned_data$LongMeth
Exp_ID <- cleaned_data$Exp_ID

#-----------------------------------
# Creating Methylation Summary Data
#-----------------------------------
GeneDepthPlant <- dcast(LongMeth,Gene*Zeroth_pos~Plant,mean, 
                        value.var = "total_RD", na.rm=TRUE)
GenePercentGroup <- create_gene_percent_x(LongPercent, 'Group', mean)
GenePercentPlant <- create_gene_percent_x(LongPercent, 'Plant', mean)

# QC
# Make sure the GenePercentX dfs are the same length as the ZoomFrame_filtered
# and in the same order
if (sum(GenePercentPlant$Zeroth_pos == ZoomFrame_filtered$Zeroth_pos) != nrow(ZoomFrame_filtered) |
    sum(GenePercentPlant$Gene == ZoomFrame_filtered$Gene) != nrow(ZoomFrame_filtered)) {
  print('Output is in a different order. Try running again.')
}


#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------
#Make version of Exp_ID file with only treated individuals
Exp_ID_Treated <- Exp_ID[Exp_ID$Group!="C",]

Output_Frame <- create_output_frame(Exp_ID_Treated, ZoomFrame_filtered,
                                    GenePercentPlant, GeneDepthPlant,
                                    GenePercentGroup, 'C')


#--------------------
# Group DMR Analysis
#--------------------
# Subset the Exp ID
Exp_ID_sub <- create_sub_ID(Exp_ID, 'Generation', 'T1')

# Run the Group DMR analysis
Output_Frame <- DMR(Output_Frame, ZoomFrame_filtered, Exp_ID, 
                    fixed = c('Group'), random = c('Plant'), 
                    colnames_of_interest, reads_threshold = 3,
                    model = 'binomial', analysis_type = 'group')

#----------------
# Individual DMR
#----------------
# Run the Individual DMR analysis
Output_Frame <- DMR(Output_Frame, ZoomFrame_filtered, Exp_ID, 
                    fixed = c('Group'), random = c('Individual'), 
                    reads_threshold = 5, control = 'C', 
                    model = 'beta-binomial', analysis_type = 'individual')

#----------------------
# Changepoint Analysis
#----------------------
# Get the potential column names to run changepoint analysis on
changepoint_cols = find_changepoint_col_options(Output_Frame)

# Create the penalties
CG_penalty <- 9; CHG_penalty <- 4; CHH_penalty <- 7



# Run the changepoint_analysis function
test1 <- changepoint_analysis(df, CG_penalty, CHG_penalty, 
                              CHH_penalty, 'Z_Zeb_Concentration_small')
# Specify the Z score column (depends on fixed effects)

#------------------------------------Write Data---------------------------------

# Write data to csv
write_csv(changepoint_output_frame, 'data/GDM_soy/changepoint_output_frame.csv')

