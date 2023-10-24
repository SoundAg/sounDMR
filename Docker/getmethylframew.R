#!/usr/bin/Rscript

library(GetoptLong)

suppressMessages(suppressWarnings(library(sounDMR)))

sample_count = 0
methyl_call_type="Dorado"
max_read_depth = 100

spec = "
This run step3 of sounDMR which creates Methylframe.

Usage: Rscript SounDMR_step2.R [options]

Options:
  <directory=s> directory path for the bedfiles.
  <sample_count=i> If you need to include the samples from a previous round, then enter the total number of samples from the previous round here. By default sample ennumeration starts with 'S'.
  <methyl_call_type=s> This determines the columsn to subset. Options include 'Megalodon', 'Dorado', 'Bonito', 'DSP''
  <max_read_depth=i> Read depth to filter out. P.S Higher read depths include regions with repetative elements. 
  <pattern=s> Regular expression pattern for the bed file for ex:'*Bedmethyl.bed'. P.S Make sure to name the bedmethyl files same as the pattern being provided. 
  <file_prefix=s> Add a prefix to files that are saved into the working directory while running the function.              
  <na_filter> Filter NAs for the samples. This would be 1/2 total number of samples by default.
  <verbose> print messages.

Contact: name@address
"
GetoptLong(spec, template_control = list(opt_width = 23))


setwd(directory)
pattern <- as.character(pattern)

methyl_beds <- list.files(path=".",pattern=pattern)

if(na_filter==0){
  na_filter = length(methyl_beds)/2
}

print(c(length(methyl_beds),na_filter, max_read_depth, sample_count))


df <- sounDMR::generate_methylframe(methyl_bed_list=methyl_beds, Sample_count = sample_count, 
Methyl_call_type=methyl_call_type, max_read_depth=max_read_depth, filter_NAs=na_filter, File_prefix=file_prefix)




