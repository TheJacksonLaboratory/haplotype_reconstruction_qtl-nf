#!/usr/bin/env Rscript
library(qtl2)
library(parallel)

################################################################################
# Perform initial calculations of genotype probabilities for genotyping error
# assessment and founder contributions with QC'd cross object
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20240118
################################################################################ 

args <- commandArgs(trailingOnly = TRUE)
cross <- args[1]
# cross <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/geno_probs/working_cross.RData"

# Load in DO cross data
# note: this cross has bad markers removed, but *not* bad or duplicate samples
load(cross)

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = working_cross, 
                          map = working_cross$pmap, 
                          error_prob = 0.002, 
                          lowmem = FALSE,
                          cores = (parallel::detectCores()/2), 
                          quiet = F)
save(pr, file = "pr_36state.RData")

# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr,
                                    cores = (parallel::detectCores()/2), 
                                    quiet = F)
save(apr, file = "pr_8state.RData")

# Estimate genotyping errors
pr_errorlod <- qtl2::calc_errorlod(cross = working_cross,
                                   probs = pr, 
                                   cores = (parallel::detectCores()/2))
save(pr_errorlod, file = "pr_errorlod.RData")

# Find best marginal genotype probability
m <- qtl2::maxmarg(probs = pr, 
                   cores = parallel::detectCores(), 
                   minprob = 0.5)
m[[1]][1:5,1:10]
m[[19]][1:5,1:10]
save(m, file = "maxmarg_numeric.RData")

# Count crossovers
nxo <- qtl2::count_xo(m, cores=parallel::detectCores())
save(nxo, file = "nxos.RData")
