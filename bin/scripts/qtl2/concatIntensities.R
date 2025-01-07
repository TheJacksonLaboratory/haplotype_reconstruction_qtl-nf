#!/usr/bin/env Rscript
library(fst)
library(parallel)
library(dplyr)

################################################################################
# Concatenate project intensity files.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241227
################################################################################ 

# testing
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/7b/0223e2bb4e3f54d4e908ac4ada1a15"
# setwd(test_dir)

# get chr X intensities objects
x_files <- list.files(pattern = "chrX")
x_int_list <- lapply(x_files, function(x) read.csv(x, skip = 3, check.names = F))
x_int_list <- lapply(x_int_list, function(z){
  # detect na intensities
  na_samples <- list()
  for(i in 1:ncol(z)){
    na_samples[i] <- all(is.na(z[,i]))
  }
  na_samples <- which(unlist(na_samples))
  
  # get rid of samples with all nas
  if(length(na_samples) == 0){
    return(z)
  } else {
    return(z[,-na_samples])
  }
})
x_int_df <- Reduce(left_join, x_int_list)

# get chr Y intensities objects
y_files <- list.files(pattern = "chrY")
y_int_list <- lapply(y_files, function(x) read.csv(x, skip = 3, check.names = F))
y_int_list <- lapply(y_int_list, function(z){
  # detect na intensities
  na_samples <- list()
  for(i in 1:ncol(z)){
    na_samples[i] <- all(is.na(z[,i]))
  }
  na_samples <- which(unlist(na_samples))
  
  # get rid of samples with all nas
  if(length(na_samples) == 0){
    return(z)
  } else {
    return(z[,-na_samples])
  }
})
y_int_df <- Reduce(left_join, y_int_list)

# .fst intensity filesq
intensity_files <- list.files(pattern = "fst")
intensity_list <- lapply(intensity_files, function(x) fst::read.fst(x))
int_df <- Reduce(left_join, intensity_list)

# save
write.csv(x_int_df, file = "chrX_intensities.csv", row.names = F, quote = F)
write.csv(y_int_df, file = "chrY_intensities.csv", row.names = F, quote = F)
fst::write.fst(int_df, path = "all_intensities.fst")
