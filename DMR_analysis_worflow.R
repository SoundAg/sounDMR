
library(sounDMR2)


#-----------------------------------------------------Part 1 : ONT data clean up and standardization---------------------------------------------------------------------------------

#-------------------------------
# Read In Gene Co-ordinates File
#-------------------------------
#It is necessary for the gene coordinates file to be in the below format to ensure the code works.  
#Chromosome | Low   | High  | Gene_name	| Strand    | Gene_length   | Adapt_Low | Adapt_High

Geneco <- read.table(file.choose(), header=TRUE, sep="\t")


#-------------------------
# Read In Methyl.bed Files
#-------------------------
#The below command works only if all the bed files that you want to work with are in the working directory.
#Make sure to change the pattern based on your methyl bed file
All_methyl_beds <- list.files(path=".",pattern="*_methyl.bed")


#-------------------------
# Create Megaframe
#-------------------------
Megalist <- generate_megaframe(methyl_bed_list=All_methyl_beds, Sample_count = 0, Methyl_call_type="",  File_prefix="")

Megaframe <- Megalist[[1]]
experimental_design_df <- Megalist[[2]]

#P.S The experimental_design starter doesn't have any information with respect to treatments, rounds etc. 
#Make sure to add it for DMR analysis

#-------------------------
# Read In Methyl.bed Files
#-------------------------
Zoomframe <- generate_zoomframe(gene_cord_df = Geneco, MFrame = Megaframe, Gene_col="Gene_name", filter_NAs=0, target_info=TRUE, gene_list = Geneco$Gene_name, File_prefix="")


#-----------------------------------------------------Part 2 : DMR analysis----------------------------------------------------------------------------------------------------------

#------------
# Clean Data
#------------

# The colnames of interest. This will be carried through the rest of the analysis
required_columns <- c('Chromosome', 'Gene', 'Position', 'Strand', 'CX', 
                          'Zeroth_pos', 'Plant')

dmr_obj <- create_dmr_obj(Zoomframe, experimental_design_df, required_columns)

#-----------------------------------
# Creating Methylation Summary Data
#-----------------------------------
GeneDepthPlant <- dcast(dmr_obj$LongMeth,Gene*Zeroth_pos~Plant,mean, 
                        value.var = "total_RD", na.rm=TRUE)
GenePercentGroup <- create_gene_percent_x(dmr_obj$LongPercent, 'Group', mean)
GenePercentPlant <- create_gene_percent_x(dmr_obj$LongPercent, 'Plant', mean)

# QC
# Make sure the GenePercentX dfs are the same length as the Zoomframe_filtered
# and in the same order
if (sum(GenePercentPlant$Zeroth_pos == dmr_obj$Zoomframe_filtered$Zeroth_pos) != nrow(dmr_obj$Zoomframe_filtered) |
    sum(GenePercentPlant$Gene == dmr_obj$Zoomframe_filtered$Gene) != nrow(dmr_obj$Zoomframe_filtered)) {
  print('Output is in a different order. Try running again.')
}


#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------

methyl_summary <- create_methyl_summary(dmr_obj, GenePercentPlant, GeneDepthPlant,
                                    GenePercentGroup, control = 'C')


#--------------------
# Group DMR Analysis
#--------------------

# Run the Group DMR analysis
methyl_summary <- DMR(methyl_summary, dmr_obj,
                    fixed = c('Group'), random = c('Plant'), 
                    required_columns, reads_threshold = 3, control = 'C',
                    model = 'binomial', analysis_type = 'group')

#----------------
# Individual DMR
#----------------
# Run the Individual DMR analysis
methyl_summary <- DMR(methyl_summary, dmr_obj, 
                    fixed = c('Group'), random = c('Individual'), 
                    required_columns, reads_threshold = 5, control = 'C', 
                    model = 'beta-binomial', analysis_type = 'individual')

#----------------------
# Changepoint Analysis
#----------------------
# Get the potential column names to run changepoint analysis on
changepoint_cols = find_changepoint_col_options(methyl_summary)

# Run the changepoint_analysis function
methyl_summary <- changepoint_analysis(methyl_summary, CG_penalty = 9, CHG_penalty = 4, 
                              CHH_penalty = 7, z_col = 'Z_Zeb_Concentration_small')



