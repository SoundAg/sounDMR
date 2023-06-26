
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
All_methyl_beds <- list.files(path=".",pattern="*methyl.bed")


#-------------------------
# Create Methylframe 
#-------------------------
#If gene_info is false in the below parameter then this returns a megaframe and if True it returns the zoomframe
Methylframe <- generate_methylframe(methyl_bed_list=All_methyl_beds, Sample_count = 0, 
                                  Methyl_call_type="Dorado", filter_NAs = 0,
                                  gene_info = FALSE, gene_coordinate_file = Geneco, Gene_column='Gene_name',
                                  target_info=FALSE, 
                                  File_prefix="Sample")


#P.S The above function creates The experimental_design starter doesn't have any information with respect to treatments, rounds etc. 
#Make sure to add it for DMR analysis



#----------------------------Part 2 : DMR analysis------------------------------

#------------
# Clean Data
#------------

experimental_design_df <- read.table(file.choose(), header=TRUE, sep=",")
dmr_obj <- create_dmr_obj(Methylframe, experimental_design_df)

#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------

methyl_summary <- create_methyl_summary(dmr_obj, control = 'C', treated = 'T',
                                        additional_summary_cols = list(c(sd, 'Group')))

# Option to subset methyl_summary
individuals_of_interest = unique(dmr_obj$experimental_design_df$Individual)
methyl_summary <- subset_methyl_summary(methyl_summary, 
                                        individuals_to_keep = individuals_of_interest)

#--------------------
# Group DMR Analysis
#--------------------

# Run the Group DMR analysis
methyl_summary_gene_info <- find_DMR(methyl_summary_gene_info, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 3, 
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
changepoint_cols = find_changepoint_col_options(methyl_summary_no_gene_info)

# The target genes of interest
target_genes <- unique(dmr_obj_no_gene_info$ZoomFrame_filtered$Gene)

# Run the changepoint_analysis function
methyl_summary_gene_info <- changepoint_analysis(methyl_summary_gene_info, CG_penalty = 9, 
                                       CHG_penalty = 4, CHH_penalty = 7, 
                                       target_genes = target_genes,
                                       save_plots = F,
                                       z_col = "Z_GroupT_small")

#----------------------
# DMR score rendering
#----------------------

DMR_score_gene_info <- sound_score(changepoint_OF = methyl_summary_gene_info, 
                         Statistic = changepoint_cols[1], 
                         Per_Change = "Treat_V_Control", CF = T,
                         other_columns=c("Control", "Estimate_GroupT_small"),
                         UserFunction = NA)

# Only run bootscore if gene info is available
DMR_boot_score <- boot_score(sound_score_obj = DMR_score_gene_info, 
                             target_gene = "AT1G01640", scoring_col_name="dmr_score2")
