#!/usr/bin/env Rscript
library(qtl2)
library(qtl2convert)
library(dplyr)
library(vroom)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
# 1) 
# Write control file for DO mice
chr <- c(1:19, "X")
write_control_file(output_file = "data/DOforqtl2.json",
                   crosstype="do",
                   description="DO_HR",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("geno", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("DO_genos/DO_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   covar_file="DO_covar.csv",
                   sex_covar="Sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="Generation", overwrite = T)

# Write control file for 4WC mice
write_control_file(output_file = "data/4WCforqtl2.json",
                   crosstype="genail4",
                   description="4WC_HR",
                   # Update lines below when these consensus genotypes actually exist; reading in this cross throws a warning
                   founder_geno_file=paste0("4WC_founder_genos/4WC_foundergeno", chr, ".csv"), 
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("4WC_genos/4WC_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   sex_covar="Sex",
                   sex_codes=list(F="Female", M="Male"), 
                   covar_file = "4WC_covar.csv",
                   crossinfo_covar = c("Generation","F","J","K","L"),
                   overwrite = T)



