#!/usr/bin/env Rscript
library(qtl2)
library(parallel)

################################################################################
# Concatenate project genotype probabilities, convert to allele probabilities, 
# calculate genotyping errors, etc.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241219
################################################################################ 

# # testing
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/a8/b534c13c57cef356e89bc79a917a9a"
# setwd(test_dir)

# get genoprobs objects
cross_file_list <- c(t(read.table("crosses.txt")))
cross_list <- lapply(cross_file_list, function(x){
  cross <- readRDS(x)
  return(cross)
})

# combine sample genotypes
new_geno <- list()
for(c in names(cross_list[[1]]$geno)){
  message(paste("Chromosome",c))
  geno <- Reduce(rbind,lapply(cross_list, function(x){
    g <- x$geno[[c]]
    return(g)
  }))
  new_geno[[c]] <- geno
}

# combine covars
new_covar <- Reduce(rbind,lapply(cross_list, function(x){
  cov <- x$covar
  return(cov)
}))

# combine crossinfos
new_crossinfo <- Reduce(rbind,lapply(cross_list, function(x){
  ci <- x$cross_info
  return(ci)
}))

# combine isfemales
new_isfemales <- unlist(lapply(cross_list, function(x) x$is_female))

# skeleton cross to replace with concatenated sample data
cross <- cross_list[[1]]
cross$covar <- new_covar
cross$cross_info <- new_crossinfo
cross$is_female <- new_isfemales

# get genoprobs objects
probs_file_list <- c(t(read.table("probs.txt")))

# read probs in
probs_list <- lapply(probs_file_list, function(x){
  pr <- readRDS(x)
  return(pr)
})

# concatenate
genoprobs <- Reduce(rbind, probs_list)

# convert to allele probs
alleleprobs <- qtl2::genoprob_to_alleleprob(probs = genoprobs, 
                                            cores = parallel::detectCores(), 
                                            quiet = F)

# save
saveRDS(cross, file = "cross.rds")
saveRDS(genoprobs, file = "genoprobs.rds", compress = "gzip")
saveRDS(alleleprobs, file = "alleleprobs.rds", compress = "gzip")