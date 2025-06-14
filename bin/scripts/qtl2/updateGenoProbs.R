#!/usr/bin/env Rscript
library(qtl2)
library(parallel)
library(dplyr)

################################################################################
# Update genotype probabilities, convert to allele probabilities, 
# calculate genotyping errors, etc.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250523
################################################################################

# # testing
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/88/c6ee71e4f7c91e94491e64f99873b8"
# setwd(test_dir)

args <- commandArgs(trailingOnly = TRUE)

covar_file <- args[1]
cross_file <- args[2]
genoprobs_file <- args[3]
alleleprobs_file <- args[4]
viterbi_file <- args[5]
kinship_file <- args[6]
marker_file <- args[7]
cross_type <- args[8]
remove_markers <- args[9]
correct_ids <- args[10]
print(args)

# Read in the files
message("Reading in files...")
covar <- read.csv(covar_file, header = TRUE, stringsAsFactors = FALSE)
if("include" %in% colnames(covar)){
    covar <- covar[covar$include == TRUE,] %>%
                dplyr::select(-include)
}
message("Covariate file read in.")
message("Reading in cross object...")
cross <- readRDS(cross_file)

# Convert logical parameters from string
remove_markers <- as.logical(remove_markers)
correct_ids <- as.logical(correct_ids)

if(remove_markers == TRUE){
    if(remove_markers == TRUE && is.null(marker_file)){
    stop("remove_markers is TRUE, but no marker file provided.")
    }
    
    message("Remove markers is TRUE, reading in bad markers...")
    bad_markers <- readRDS(marker_file)

    message("Dropping bad markers...")
    cross <- qtl2::drop_markers(cross, bad_markers)

    message("Estimating new genotype probabilities...")
    # Calculate genotype probs
    genoprobs <- qtl2::calc_genoprob(cross = cross, 
                            map = cross$pmap, 
                            error_prob = 0.002,
                            lowmem = FALSE,
                            cores = (parallel::detectCores()/2), 
                            quiet = F)
    
    message("Converting to allele probabilities...")
    # convert to allele probs
    alleleprobs <- qtl2::genoprob_to_alleleprob(probs = genoprobs, 
                                                cores = parallel::detectCores(), 
                                                quiet = F)

    # calculate kinship matrix
    k <- calc_kinship(genoprobs, type = "loco", quiet = FALSE, cores = 0)

    # calculate viterbi
    m <- qtl2::maxmarg(probs = genoprobs, minprob=0.5, map = cross$pmap, quiet = T)


} else{
    message("Remove markers is FALSE, reading probability files...")
    # Read in files that we need if correct_ids is TRUE and remove_markers is FALSE
    genoprobs <- readRDS(genoprobs_file)
    alleleprobs <- readRDS(alleleprobs_file)
    m <- readRDS(viterbi_file)
    k <- readRDS(kinship_file)
}

if(correct_ids == TRUE){
        message("Correcting IDs...")
        if("correct_id" %in% colnames(covar) == FALSE){
            stop("correct_ids is TRUE, but no correct_id column in covar file.")
        }
        # Make correct_ids vector
        correct_ids <- covar$correct_id
        names(correct_ids) <- covar$id

        # Correct IDs
        cross <- qtl2::replace_ids(cross, correct_ids)
        genoprobs <- qtl2::replace_ids(genoprobs, correct_ids)
        alleleprobs <- qtl2::replace_ids(alleleprobs, correct_ids)
        m <- qtl2::replace_ids(m, correct_ids)
        k <- qtl2::replace_ids(k, correct_ids)

        message("IDs corrected.")
        message("Saving updated cross, genotype probabilities, allele probabilities, and viterbi files...")
        saveRDS(cross, file = "cross.rds")
        saveRDS(genoprobs, file = "genoprobs.rds", compress = TRUE)
        saveRDS(alleleprobs, file = "alleleprobs.rds", compress = TRUE)
        saveRDS(m, file = "maxmarg.rds", compress = TRUE)
        saveRDS(k, file = "kinship.rds", compress = TRUE)
        write.csv(covar, file = "covar.csv", row.names = FALSE, quote = FALSE)


} else {
        message("rerun = FALSE, remove_markers = FALSE, and correct_ids = FALSE. No changes required to haplotype reconstructions.")
        message("Saving updated cross, genotype probabilities, allele probabilities, and viterbi files...")
        saveRDS(cross, file = "cross.rds")
        saveRDS(genoprobs, file = "genoprobs.rds", compress = TRUE)
        saveRDS(alleleprobs, file = "alleleprobs.rds", compress = TRUE)
        saveRDS(m, file = "maxmarg.rds", compress = TRUE)
        saveRDS(k, file = "kinship.rds", compress = TRUE)
        write.csv(covar, file = "covar.csv", row.names = FALSE, quote = FALSE)
}








