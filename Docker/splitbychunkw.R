#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


library(sounDMR)
result <- sounDMR::split_by_chunk(args[1],args[2],args[3])
print(result)
