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

Usage: Rscript step4_sounDMRw.R [options]

Options:
  <directory=s> Directory path for the methylframe/megaframe. Hint:chunks
  <Megaframe_pattern=s> Pattern of megaframe file.
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

print(directory)
print(exp_design)

experimental_design_df <- fread(exp_design, header = TRUE, sep=",")

base_dir <- directory
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
process_subdir <- function(subdir) {
  mf <- paste(subdir, "megaframe", sep="/")
  # mf: /ftmp/batch1_4/chr_chr1/chunks/chunk1_0_250000/megaframe
  groups <- strsplit(subdir, split = "/")
  chrmall <- groups[[1]][6]
  # chrmall: chunk1_0_250000
  chk <- groups[[1]][4]
  # chk: chr_chr1
  chrm <- gsub("chr_", "", chk)
  file_pre <- paste(file_prefix, chrm, chrmall, sep = "_")
  # file_pre: FP_chr1_chunk1_0_250000
  setwd(mf)
  if (length(Sys.glob(Megaframe_pattern, dirmark = FALSE))>1) {
  stop("There is more than one pattern matching for Megaframe_pattern")
  }
  Methylframe <- fread(Sys.glob(Megaframe_pattern, dirmark = FALSE)[1], header=TRUE, sep=",")
  #create dmrobj
  dmr_obj <- sounDMR::create_dmr_obj(Methylframe, experimental_design_df)
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


