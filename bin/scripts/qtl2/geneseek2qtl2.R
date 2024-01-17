#!/usr/bin/env Rscript

# convert GeneSeek FinalReport files to format for R/qtl2
#
# - creates one genotype CSV file for each chromosome
#
# - also creates 2 files containing the two channels of SNP intensities for markers on the X and Y chr
#   (these are useful for verifying the sex of the mice)

# load required packages
library(data.table)
library(qtl2convert)
library(vroom)
library(parallel)
library(dplyr)
library(magrittr)
library(fst)
# library(microbenchmark)
################################################################################
# Parse Neogen FinalReport files and create intensity files and genotypes.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20240102
################################################################################ 
args <- commandArgs(trailingOnly = TRUE)

# test dir
# test_dir <- "/fastscratch/QC_HAP_outputDir/work/8b/61d065d83b4d93a5c58eb3579e56a5"

# file containing allele codes for GigaMUGA data
# codefile <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data/GM_allelecodes.csv"
codefile <- args[1]

# input files with GigaMUGA genotypes
# ifiles <- file.path(test_dir,"finalreportlist.txt")
ifiles <- args[2]

# metadata
metadata <- args[3]
# metadata <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/covar_files/DO_ESC_covar.csv"
meta <- read.csv(metadata)

# read genotype codes
codes <- data.table::fread(codefile, skip = 3)

# replacement function for qtl2convert::cbind_smother
# instead of condensing repeated samples:
# records the name of the duplicate (original sample name + a three digit randomly generated "dup" flag) and appends the repeated genotypes
combineGenos <- function(big_geno, small_geno){
  if(any(colnames(small_geno) %in% colnames(full_geno))){
    rep_samples <- unlist(lapply(colnames(small_geno), function(x){
      flag <- paste0(runif(3,0,9) %/% 1, collapse = "")
      paste0(x,"_dup_",flag)
    }))
    update_meta <- data.frame(colnames(small_geno), rep_samples)
    colnames(update_meta) <- c("original","duplicated")
    colnames(small_geno) <- rep_samples
  } else {
    update_meta <- data.frame(NA,NA)
    colnames(update_meta) <- c("original","duplicated")
  }
  new_big_geno <- cbind(big_geno, small_geno)
  return(list(new_big_geno, update_meta))
}

# read in the FinalReport file table from nextflow
finalreports <- read.table(ifiles)

# create vector of FinalReports from the table
ifiles <- c(as.matrix(t(finalreports)))
# ifiles <- paste(test_dir,c(as.matrix(t(finalreports))), sep = "/")

# initialize full genotype matrix
full_geno <- NULL

# initialize sex chromosome intensity matrices
cXint <- cYint <- NULL

# initialize table of duplicated samples across FinalReports
duplicate_table <- NULL

# all marker intensities
allints <- NULL

# pull out genotypes and sex chromosome intensities from FinalReports
for(ifile in ifiles) {
    cat(" -File:", ifile, "\n")
    cat(" -Reading data\n")
    
    # read in the FinalReport compressed file
    g <- vroom::vroom(file = ifile, skip = 9,
                      # n_max = 100000,
                      num_threads = parallel::detectCores()/2, progress = F)
    
    # pull the sample names from the qtl2 covar file
    meta_sample_names <- unique(meta$id)
    
    # pull the sample names from the FinalReport file
    sample_names <- unique(g$`Sample ID`)
    # print(sample_names)
    # match the two sources
    matches <- sample_names[which(sample_names %in% meta_sample_names)]
    print(matches)


    # are *any* samples from metadata are in the finalreport file?
    if(length(matches) == 0){
      next
    }

    # are *all* samples from metadata in this finalreport file?
    if(length(matches) < length(meta_sample_names)){
      cat(paste0(" NOTE: Some samples in supplied metadata not present in", ifile,"\n"))
      cat(paste0(" Filtering", ifile," to matching samples\n"))
      g <- g[which(g$`Sample ID` %in% matches),]
      g <- g[which(g$`SNP Name` %in% codes$marker),]

    # all samples from metadata in the finalreport file
    } else {
      cat(paste0(" NOTE: All samples in supplied metadata present in", ifile,"\n"))
      g <- g[which(g$`SNP Name` %in% codes$marker),]
    }

    # create a matrix to contain the genotypes
    samples <- unique(g$`Sample ID`)
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes$marker, samples)

    # fill in the matrix
    cat(" -Reorganizing data\n")
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (g$`Sample ID`==samples[i])
        geno[g[wh,]$`SNP Name`,i] <- paste0(g[wh,]$`Allele1 - Forward`,
                                            g[wh,]$`Allele2 - Forward`)
    }

    # Encode sample genotypes using the allele codes
    cat(" -Encode genotypes\n")
    geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))

    # is this the first/only finalreport file?
    if(is.null(full_geno)) {
        full_geno <- geno

    # if not, need to combine genotypes from multiple finalreport files
    } else {
      # if there are any samples in both genotype matrices, amend the name
      combineGeno_out <- combineGenos(full_geno, geno)
      full_geno <- combineGeno_out[[1]]
      dup_table <- combineGeno_out[[2]]

      # were there any duplicates?
      if(nrow(dup_table[complete.cases(dup_table),]) != 0){
        if(is.null(duplicate_table)){
          # write duplicates to empty object
          duplicate_table <- dup_table
          
        } else {
          # append the duplicates to existing duplicates
          duplicate_table <- rbind(duplicate_table,dup_table)
        }
        
        # replace the names in the genotype file
        g <- g %>%
          dplyr::rename(original = `Sample ID`) %>%
          dplyr::left_join(., dup_table) %>%
          dplyr::rename(`Sample ID` = duplicated) %>%
          dplyr::select(`SNP Name`,`Sample ID`, everything(), -original)
        samples <- unique(g$`Sample ID`)
      }
    }
    
    # grab X and Y intensities
    cat(" -Grab X and Y intensities\n")
    gX <- g[which(g$`SNP Name` %in% codes[codes$chr=="X"]$marker),]
    gY <- g[which(g$`SNP Name` %in% codes[codes$chr=="Y"]$marker),]

    cX <- matrix(nrow=sum(codes$chr=="X"),
                 ncol=length(samples))
    dimnames(cX) <- list(codes[which(codes$chr == "X"),]$marker, samples)

    cY <- matrix(nrow=sum(codes$chr=="Y"),
                 ncol=length(samples))
    dimnames(cY) <- list(codes[which(codes$chr == "Y"),]$marker, samples)

    for(i in seq(along=samples)) {
      if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
      wh <- (gX[,"Sample ID"]==samples[i])
      cX[gX[wh,]$`SNP Name`,i] <- (gX$X[wh] + gX$Y[wh])/2

      wh <- (gY[,"Sample ID"]==samples[i])
      cY[gY[wh,]$`SNP Name`,i] <- (gY$X[wh] + gY$Y[wh])/2
      # cY[gY[wh,"SNP Name"],i] <- (gY$X[wh] + gY$Y[wh])/2
    }

    if(is.null(cXint)) {
      cXint <- cX
      cYint <- cY
    } else {
      # if any columns in both, use those from second set
      cXint <- cbind(cXint, cX)
      cYint <- cbind(cYint, cY)
    }
    
    # pull the intensities and new sample names and write to intensity object
    int <- g %>%
      dplyr::select(`SNP Name`,`Sample ID`,X,Y) %>%
      dplyr::rename(snp = `SNP Name`,
                    sample = `Sample ID`) %>%
      tidyr::pivot_longer(-c(snp,sample),names_to = "channel") %>%
      tidyr::pivot_wider(names_from = sample, values_from = value) %>%
      data.frame()
    
    if(is.null(allints)){
      allints <- int
    } else {
      allints <- dplyr::left_join(allints, int)
    }
    
    
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(cbind(marker=rownames(cXint), cXint),
                       "chrXint.csv",
                       "X chr intensities",
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint), cYint),
                       "chrYint.csv",
                       "Y chr intensities",
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chrom in c(1:19,"X","Y","M")) {
    # mar <- codes[codes$chr==chr,"marker"]
    mar <- codes %>%
      dplyr::filter(chr == chrom) %>%
      dplyr::select(marker)
    g <- full_geno[rownames(full_geno) %in% mar$marker,]
    qtl2convert::write2csv(cbind(marker=rownames(g), g),
                           paste0("geno", chrom, ".csv"),
                           paste0("genotypes for chr ", chrom),
                           overwrite=TRUE)
}

# write to fst file, maximally compressed
write_fst(allints, "intensities.fst", compress=100)

# write new metadata file
if(!is.null(duplicate_table)){
  duplicate_covar <- duplicate_table %>%
    dplyr::rename(id = original) %>%
    dplyr::left_join(., meta) %>%
    dplyr::select(-id) %>%
    dplyr::rename(id = duplicated) %>%
    dplyr::bind_rows(meta,.)
  write.csv(duplicate_covar, file = "covar.csv", quote = F, row.names = F)
} else {
  write.csv(meta, file = "covar.csv", quote = F, row.names = F)
}








