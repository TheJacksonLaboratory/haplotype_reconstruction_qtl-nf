#!/usr/bin/env Rscript
library(qtl2)
library(fst)
################################################################################
# Perform marker and sample QC using cross object and intensities upstream in 
# haplotype reconstruction pipeline.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20240112
################################################################################ 
# test_dir <- "/fastscratch/QC_HAP_outputDir/work/b8/fdfd8ce7bdd497445041b0ac27a8b7"
# setwd(test_dir)
args <- commandArgs(trailingOnly = TRUE)

# import cross object
cross <- args[1]
# cross <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/geno_probs/preQC_cross.RData"
load(cross)

# import intensities
intensities <- args[2]
# intensities <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/qtl2genos/intensities.fst"


# Reordering genotypes so that most common allele in founders is first
for(chr in seq_along(cross$founder_geno)) {
  fg <- cross$founder_geno[[chr]]
  g <- cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  cross$founder_geno[[chr]] <- fg
  cross$geno[[chr]] <- g
}

# percent missing markers
percent_missing_marker <- qtl2::n_missing(cross, "marker", "prop")*100

# markers missing in >10% of samples
missing_markers <- percent_missing_marker[percent_missing_marker > 10]

# remove those markers
working_cross <- qtl2::drop_markers(cross, markers = names(missing_markers))

# percent missing genotypes at sample level
# original cross
percent_missing_ind_cross <- qtl2::n_missing(cross, "ind", "prop")*100

# Sample Duplicates
cg <- compare_geno(working_cross, cores=0)

# Sex checks
# Reading in all probe intensities
int <- fst::read.fst(intensities)
int <- int[int$snp %in% qtl2::marker_names(working_cross),]
int <- int[seq(1, nrow(int), by=2),-(1:2)] + int[-seq(1, nrow(int), by=2),-(1:2)]
int <- int[,which(colnames(int) %in% qtl2::ind_ids(cross))]

# # sample genotypes
# g <- do.call("cbind", cross$geno[1:19])
# 
# # founder genotypes
# fg <- do.call("cbind", cross$founder_geno[1:19])
# 
# # find sample genotypes with markers with no founder genotype?
# g <- g[,colSums(fg==0)==0]
# 
# # find markers with missing genotypes in the founders
# fg <- fg[,colSums(fg==0)==0]
# fgn <- colSums(fg==3)

# save objects to send to markdown
save(cg, # sample duplicates
     int, # marker intensities
     percent_missing_ind_cross, # percent of missing genotypes per sample
     percent_missing_marker, # percent of samples missing each marker genotype; 
     # *these markers are removed in working_cross and intensities (int)*
     working_cross,
     file = "QC_1.RData")

# save just the cross object for the genoprobs process
save(working_cross, file = "working_cross.RData")