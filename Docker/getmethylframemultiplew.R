#!/usr/bin/Rscript

library(GetoptLong)

suppressMessages(suppressWarnings(library(sounDMR)))

sample_count = 0
methyl_call_type="Dorado"
max_read_depth = 100

spec = "
This run step3 of sounDMR which creates Methylframe. 
This will be applied over all directories inside given directory.

Usage: Rscript SounDMR_step3.R [options]

Options:
  <directory=s> directory path for the chunks.
  <sample_count=i> If you need to include the samples from a previous round, then enter the total number of samples from the previous round here. By default sample ennumeration starts with 'S'.
  <methyl_call_type=s> This determines the columsn to subset. Options include 'Megalodon', 'Dorado', 'Bonito', 'DSP''
  <max_read_depth=i> Read depth to filter out. P.S Higher read depths include regions with repetative elements. 
  <file_prefix=s> Add a prefix to files that are saved into the working directory while running the function.              
  <na_filter> Filter NAs for the samples. This would be 1/2 total number of samples by default.
  <verbose> print messages.

Contact: name@address
"
GetoptLong(spec, template_control = list(opt_width = 23))

base_dir <- directory
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
process_subdir <- function(subdir) {
  print(subdir)
  setwd(subdir)
  methyl_beds <- list.files(path=".", pattern='*.bed')
  if(na_filter==0){
  na_filter = length(methyl_beds)/2
  }
  df <- sounDMR::generate_methylframe(methyl_bed_list=methyl_beds, Sample_count = sample_count, Methyl_call_type=methyl_call_type, max_read_depth=max_read_depth, filter_NAs=na_filter, File_prefix=file_prefix)
}

lapply(subdirs, process_subdir)






