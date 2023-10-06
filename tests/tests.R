library(testthat)
library(sounDMR)


test_that("split_by_chromosome works correctly", {
  # read test org
  good_bed_dfs <- list()
  for (i in 0:5) {
    fn <- paste0("syntheticorg.", i,"bed.")        
    idx <- paste0("chr", i)
    good_bed_dfs[[idx]] <- read.table(fn, sep = "\t", header = FALSE)
  }
  # run function
  fn <- paste0(getwd(), "syntheticorg.bed")
  result <- sounDMR::split_by_chromosome(fn)
  # get results
  result_bed_dfs <- list()
  for (i in 0:5) {
    fn <- paste0("chr_SL4.0ch0", i, "/syntheticorg.bed")    
    idx <- paste0("chr", i)
    result_bed_dfs[[idx]] <- read.table(fn, sep = "\t", header = FALSE)
    # clean the directory
    unlink(paste0("chr_SL4.0ch0", i), recursive = TRUE)
  }
  # compare
  for (i in 0:5) {
    idx <- paste0("chr", i)
    expect_equal(good_bed_dfs[idx], result_bed_dfs[idx])
  }
})
