
library(sounDMR2)

#-----------------Part 1 : ONT data clean up and standardization----------------

#-------------------------------
# Read In Gene Co-ordinates File
#-------------------------------
#It is necessary for the gene coordinates file to be in the below format to ensure the code works.  
#Chromosome | Low   | High  | Gene_name	| Strand    | Gene_length   | Adapt_Low | Adapt_High

Geneco <- read.table(file.choose(), header=TRUE, sep=",")


#-------------------------
# Read In Methyl.bed Files
#-------------------------
#The below command works only if all the bed files that you want to work with are in the working directory.
#Make sure to change the pattern based on your methyl bed file
All_methyl_beds <- list.files(path=".",pattern="*_methyl.bed")


#-------------------------
# Create Megaframe
#-------------------------
Megaframe <- generate_megaframe(methyl_bed_list=All_methyl_beds, Sample_count = 0, 
                                Methyl_call_type="DSP",  File_prefix="Sample")


#P.S The above function creates The experimental_design starter doesn't have any information with respect to treatments, rounds etc. 
#Make sure to add it for DMR analysis


#-------------------------
# Create Zoomframe
#-------------------------
Zoomframe <- generate_zoomframe(gene_cord_df = Geneco, MFrame = Megaframe, 
                                Gene_col="Gene_name", filter_NAs=2, 
                                target_info=TRUE, gene_list = Geneco$Gene_name, 
                                File_prefix="Sample")


#----------------------------Part 2 : DMR analysis------------------------------

#------------
# Clean Data
#------------

dmr_obj <- create_dmr_obj(Zoomframe, experimental_design_df)

#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------

methyl_summary <- create_methyl_summary(dmr_obj, control = 'C')

# Option to subset methyl_summary
indiduals_of_interest = c()
methyl_summary_subset <- subset_methyl_summary(methyl_summary, 
                                               individuals_to_keep = individuals_of_interest)

#--------------------
# Group DMR Analysis
#--------------------

# Run the Group DMR analysis
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Plant'), reads_threshold = 3, 
                           control = 'C', model = 'binomial', 
                           analysis_type = 'group')

#----------------
# Individual DMR
#----------------
# Run the Individual DMR analysis
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 5, 
                           control = 'C', model = 'beta-binomial', 
                           analysis_type = 'individual')

#----------------------
# Changepoint Analysis
#----------------------
# Get the potential column names to run changepoint analysis on
changepoint_cols = find_changepoint_col_options(methyl_summary)

# The target genes of interest
target_genes <- c()

# Run the changepoint_analysis function
methyl_summary <- changepoint_analysis(methyl_summary, CG_penalty = 9, 
                                       CHG_penalty = 4, CHH_penalty = 7, 
                                       target_genes = target_genes,
                                       save_plots = F,
                                       z_col = "Z_GroupT_small")

#----------------------
# DMR score rendering
#----------------------

DMR_score <- sound_score(changepoint_OF = methyl_summary, 
                         Statistic = changepoint_cols[1], 
                         Per_Change = "Treat_V_Control", 
                         other_columns=c("Control", "Estimate_GroupT_small"))


boot_score(sound_score_obj = DMR_score, target_gene = "AT1G01640", scoring_col_name="dmr_score2")
