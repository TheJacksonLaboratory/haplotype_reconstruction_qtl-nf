#!/usr/bin/env Rscript

################################################################################
# Make a .csv manifest for haplotype reconstruction from MUGA-series files.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241112
################################################################################

library(tidyverse)

# hr-nf directory
hr_nf_dir <- '/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf'

# output directory
out_dir <- file.path(hr_nf_dir,"sample_sheets")

# neogen files
finalreport_file <- list.files(hr_nf_dir, pattern = "FinalReport", recursive = T, full.names = T)[1:8]

# project ids
project_id <- unlist(lapply(finalreport_files, function(x) strsplit(x,"/")[[1]][[9]]))

# covar files
covar_file <- unlist(lapply(project_ids, function(x){
  list.files(file.path(hr_nf_dir, "projects", x, "covar_files"), recursive = T, full.names = T)
}))

# cross type
cross_type <- rep("do", length(covar_file))

# make .csv manifest
manifest <- data.frame(finalreport_file, project_id, covar_file, cross_type)

# write manifest
write.csv(manifest, file = file.path(out_dir,"test_sheet.csv"), quote = F, row.names = F)

