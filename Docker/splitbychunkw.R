#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)


library(sounDMR)
chuncksize <- as.integer(args[3])
result <- sounDMR::split_by_chunk(args[1], args[2], chuncksize)
print(result)
