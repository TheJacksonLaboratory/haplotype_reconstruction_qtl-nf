#!/usr/bin/env Rscript
library(qtl2)
library(fst)
################################################################################
# Perform marker and sample QC using cross object and intensities upstream in 
# haplotype reconstruction pipeline.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241113
################################################################################ 

args <- commandArgs(trailingOnly = TRUE)

## testing
# cross <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/results/DO_ESC/qtl2genos/preQC_cross.RData"
# intensities <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/results/DO_ESC/intensities/intensities.fst"

# cross object
cross <- args[1]

# intensities in .fst file
intensities <- args[2]

# load cross object
load(cross)

# percent missing markers
percent_missing_marker <- qtl2::n_missing(cross, "marker", "prop")*100

# markers missing in >10% of samples
missing_markers <- percent_missing_marker[percent_missing_marker > 10]

# percent missing genotypes at sample level
# original cross
percent_missing_ind_cross <- qtl2::n_missing(cross, "ind", "prop")*100
inds_with_missing_markers <- which(percent_missing_ind_cross > 10)

# Sample Duplicates
cg <- compare_geno(cross, cores=0)

# Reading in all probe intensities and filtering to appropriate markers and individuals
int <- fst::read.fst(intensities)
int <- int[int$snp %in% qtl2::marker_names(cross),]
int <- int[seq(1, nrow(int), by=2),-(1:2)] + int[-seq(1, nrow(int), by=2),-(1:2)]
int <- int[,which(colnames(int) %in% qtl2::ind_ids(cross))]

# remove markers with > 10% missingness
working_cross <- qtl2::drop_markers(cross, markers = names(missing_markers))
if(length(inds_with_missing_markers) > 0){
  working_cross <- subset(working_cross, qtl2::ind_ids(cross)[!qtl2::ind_ids(cross) %in% inds_with_missing_markers])
}

# sample genotypes
g <- do.call("cbind", cross$geno[1:19])

# founder genotypes
fg <- do.call("cbind", cross$founder_geno[1:19])

# find markers with markers with no founder genotype?
null_founder_genos <- any(colSums(fg) == 0)

# rename cross for clarity
original_cross <- cross

if(null_founder_genos){
  null_founder_markers <- colnames(fg[,which(colSums(fg) == 0)])
  working_cross <- qtl2::drop_markers(cross, markers = null_founder_markers)
  
  # save objects to send to markdown
  save(cg, # sample duplicates
       int, # marker intensities
       percent_missing_ind_cross, # percent of missing genotypes per sample
       percent_missing_marker, # percent of samples missing each marker genotype; 
       null_founder_markers, # any markers where there aren't founder genotypes
       working_cross,
       original_cross,
       file = "QC_1.RData")
  
  # save just the cross object with putative good markers and samples for the genoprobs process
  save(working_cross, file = "working_cross.RData")
  
} else {
  
  # save objects to send to markdown
  save(cg, # sample duplicates
       int, # marker intensities
       percent_missing_ind_cross, # percent of missing genotypes per sample
       percent_missing_marker, # percent of samples missing each marker genotype; 
       working_cross,
       original_cross,
       file = "QC_1.RData")
  
  # save just the cross object with putative good markers and samples for the genoprobs process
  save(working_cross, file = "working_cross.RData")
  
}

