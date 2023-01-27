#!/usr/bin/env

#-------------------------------------------------------------------------------

# Name: Optimized_DMR
# Author: Tom Cairns
# Purpose: This script is the workflow for the DMR analysis.

#-------------------------------------------------------------------------------

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

dmr_obj <- create_dmr_obj(ZoomFrame, Exp_ID, colnames_of_interest)

#-----------------------------------
# Creating Methylation Summary Data
#-----------------------------------
GeneDepthPlant <- dcast(dmr_obj$LongMeth,Gene*Zeroth_pos~Plant,mean, 
                        value.var = "total_RD", na.rm=TRUE)
GenePercentGroup <- create_gene_percent_x(dmr_obj$LongPercent, 'Group', mean)
GenePercentPlant <- create_gene_percent_x(dmr_obj$LongPercent, 'Plant', mean)

# QC
# Make sure the GenePercentX dfs are the same length as the ZoomFrame_filtered
# and in the same order
if (sum(GenePercentPlant$Zeroth_pos == dmr_obj$ZoomFrame_filtered$Zeroth_pos) != nrow(dmr_obj$ZoomFrame_filtered) |
    sum(GenePercentPlant$Gene == dmr_obj$ZoomFrame_filtered$Gene) != nrow(dmr_obj$ZoomFrame_filtered)) {
  print('Output is in a different order. Try running again.')
}


#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------

Output_Frame <- create_output_frame(dmr_obj, GenePercentPlant, GeneDepthPlant,
                                    GenePercentGroup, 'C')


#--------------------
# Group DMR Analysis
#--------------------

# Run the Group DMR analysis
Output_Frame <- DMR(Output_Frame, dmr_obj
                    fixed = c('Group'), random = c('Plant'), 
                    colnames_of_interest, reads_threshold = 3,
                    model = 'binomial', analysis_type = 'group')

#----------------
# Individual DMR
#----------------
# Run the Individual DMR analysis
Output_Frame <- DMR(Output_Frame, dmr_obj, 
                    fixed = c('Group'), random = c('Individual'), 
                    reads_threshold = 5, control = 'C', 
                    model = 'beta-binomial', analysis_type = 'individual')

#----------------------
# Changepoint Analysis
#----------------------
# Get the potential column names to run changepoint analysis on
changepoint_cols = find_changepoint_col_options(Output_Frame)

# Run the changepoint_analysis function
Output_Frame <- changepoint_analysis(Output_Frame, CG_penalty = 9, CHG_penalty = 4, 
                              CHH_penalty = 7, z_col = 'Z_Zeb_Concentration_small')

#------------------------------------Write Data---------------------------------

# Write data to csv
write_csv(changepoint_output_frame, 'data/changepoint_output_frame.csv')

