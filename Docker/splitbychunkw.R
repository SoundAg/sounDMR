#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


library(sounDMR)
chuncksize <- as.integer(args[2])
result <- sounDMR::split_by_chunk(args[1], chuncksize, args[3])
print(result)
