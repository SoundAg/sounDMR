#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(sounDMR)
result <- sounDMR::split_by_chromosome(args[1])
print(result)
