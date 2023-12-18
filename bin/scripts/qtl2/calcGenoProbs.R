#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(ggplot2)
library(parallel)
library(ggbeeswarm)
library(future)
library(furrr)


################################################################################
# Perform initial calculations of genotype probabilities for genotyping error
# assessment and founder contributions
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20231215
################################################################################ 

args <- commandArgs(trailingOnly = TRUE)
cross <- args[1]

# Load in DO cross data
load(cross)

# Drop null markers
cross <- qtl2::drop_nullmarkers(cross)

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2), quiet = F)
save(pr, file = "pr_36state.RData")

# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr,
                                    cores = (parallel::detectCores()/2), quiet = F)
save(apr, file = "pr_8state.RData")

# Estimate genotyping errors
pr_errorlod <- qtl2::calc_errorlod(cross = cross,
                                   probs = pr, 
                                   cores = (parallel::detectCores()/2))
save(pr_errorlod, file = "pr_errorlod.RData")

# Find best marginal genotype probability
m <- qtl2::maxmarg(probs = pr, 
                   cores = parallel::detectCores(), 
                   minprob = 0.5)
save(m, file = "maxmarg_numeric.RData")

# Count crossovers
nxo <- qtl2::count_xo(m, cores=parallel::detectCores())
save(nxo, file = "nxos.RData")
