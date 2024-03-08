#!/usr/bin/Rscript

suppressMessages(library(GetoptLong))
suppressMessages(suppressWarnings(library(sounDMR)))

control_value = 'C' 
treated_value = 'T'
additional_sum_col = 'Group'
stat = 'sd'

spec = "
This creates a dmr object and methyl syummary. the methyl summary will be saved in the directory.  
This will be applied over all directories inside given directory.

Usage: Rscript SounDMR_step3.R [options]

Options:
  <directory=s> Directory path for the methylframe/megaframe. Hint:chunks
  <Megaframe_suffix=s> Suffix of megaframe file.
  <exp_design=s> Input experimental design file. 
  <control_value=s> Default is C. This takes in the variable for control in the additional summary column. For example if your Group is T vs C, then control column would be C and treated column will be T, additional_sum_col will be 'Group'
  <treated_value=s> Default is T. This takes in the variable for control in the additional summary column.
  <additional_sum_col=s> Default is Group. It could be Treatment, Group, Treatment type etc. based on the desig, change the column values accordingly
  <stat=s>  Default is sd (standard deviation). Apart from mean, enter other types of statistic to obtain from the groups
  <file_prefix=s> Add a prefix to files that are saved into the working directory while running the function.
  <verbose> print messages.

Contact: name@address
"
GetoptLong(spec, template_control = list(opt_width = 23))

experimental_design_df <- fread(exp_design, header = TRUE, sep=",")

base_dir <- directory
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
process_subdir <- function(subdir) {
  print(subdir)
  
  groups <- strsplit(subdir, split = "/")
  chrmall <- groups[[1]][6]
  chk <- groups[[1]][4]
  chrm <- gsub("chr_", "", chk)
  file_pre <- paste(file_prefix, chrm, chrmall, sep = "_")
  
  setwd(subdir)
  
  #Sebastian_to fix
  #File to choose for test: /efs/dmr/batch1-4/chr_chr1/chunks/chunk1_0_250000/megaframe/chunk1_0_250000_chr1_B1_4_combined_MegaFrame.csv
  Methylframe <- fread(paste('*',Megaframe_suffix), header=TRUE, sep=",")
  #create dmrobj
  dmr_obj <- soundMR::create_dmr_obj(Methylframe, experimental_design_df)
  #create methyl summary
  methyl_summary <- soundMR::create_methyl_summary(dmr_obj, 
                                          control = control_value, 
                                          treated = treated_value,
                                          additional_summary_cols = list(
                                            c(stat, additional_sum_col)))
  
  
  write.csv(methyl_summary,paste(file_pre,"methyl_sum.csv", sep="_"), row.names = F)
  rm(dmrobj,methyl_summary)
}
lapply(subdirs, process_subdir)


