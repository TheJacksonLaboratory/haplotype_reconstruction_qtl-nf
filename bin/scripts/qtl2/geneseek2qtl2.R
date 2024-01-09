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
test_dir <- "/fastscratch/QC_HAP_outputDir/work/7f/659003a192d7f73ed3590964d5fe96"

# file containing allele codes for GigaMUGA data
codefile <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data/GM_allelecodes.csv"
# codefile <- args[1]

# input files with GigaMUGA genotypes
ifiles <- file.path(test_dir,"finalreportlist.txt")
# ifiles <- args[2]

# metadata
# metadata <- args[3]
metadata <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/covar_files/DO_ESC_covar.csv"
meta <- read.csv(metadata)

# read genotype codes
codes <- data.table::fread(codefile, skip = 3)

full_geno <- NULL
cXint <- cYint <- NULL
duplicate_table <- NULL

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


finalreports <- read.table(ifiles)
# ifiles <- c(as.matrix(t(finalreports)))
ifiles <- paste(test_dir,c(as.matrix(t(finalreports))), sep = "/")

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
    
    # match the two sources
    matches <- sample_names[which(sample_names %in% meta_sample_names)]
  
    
    # # find the samples from covariate file in the FinalReport file
    # name_check <- Reduce(rbind, 
    #        lapply(meta_sample_names, 
    #               function(x) data.frame(x, length(grep(x = sample_names, 
    #                                                     pattern = x)))))
    # colnames(name_check) <- c("sample","matches")
    
    # test if *some* samples from metadata are in the finalreport file
    if(length(matches) == 0){
      next
    }
    
    if(length(matches) < length(meta_sample_names)){
      cat(paste0(" NOTE: Some samples in supplied metadata not present in", ifile,"\n"))
      cat(paste0(" Filtering", ifile," to matching samples\n"))
      cat(paste0("Matching samples:\n"))
      print(matches)
      g %<>%
        dplyr::filter(`SNP Name` %in% codes$marker,
                      `Sample ID` %in% matches)

    } else {
      cat(paste0(" NOTE: All samples in supplied metadata present in", ifile,"\n"))
      g %<>%
        dplyr::filter(`SNP Name` %in% codes$marker)
    }
    
    # NOTE: may need to revise the IDs in the 2nd column
    samples <- unique(g$`Sample ID`)

    # matrix to contain the genotypes
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes$marker, samples)

    # fill in matrix
    cat(" -Reorganizing data\n")
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (g$`Sample ID`==samples[i])
        geno[g[wh,]$`SNP Name`,i] <- paste0(g[wh,]$`Allele1 - Forward`,
                                            g[wh,]$`Allele2 - Forward`)
    }

    cat(" -Encode genotypes\n")
    geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))

    if(is.null(full_geno)) {
        full_geno <- geno
    } else {
        # if there are any samples in both genotype matrices, amend the name and
        # update covar file
        combineGeno_out <- combineGenos(full_geno, geno)
        full_geno <- combineGeno_out[[1]]
        if(is.null(duplicate_table)){
          duplicate_table <- combineGeno_out[[2]]
        } else {
          duplicate_table <- rbind(duplicate_table,combineGeno_out[[2]])
          g <- dplyr::left_join(g, duplicate_table %>%
                                  dplyr::rename(`Sample ID` = original)) %>%
            dplyr::select(-c(`Sample ID`)) %>%
            dplyr::mutate(`Sample ID` = duplicated) %>%
            dplyr::select(`SNP Name`,`Sample ID`, everything())
        }
    }
    
    # grab X and Y intensities
    cat(" -Grab X and Y intensities\n")
    # gX <- g[g[,"SNP Name"] %in% codes[codes$chr=="X","marker"],]
    gX <- g %>%
        dplyr::filter(`SNP Name` %in% codes[which(codes$chr == "X"),]$marker)
    # gY <- g[g[,"SNP Name"] %in% codes[codes$chr=="Y","marker"],]
    gY <- g %>%
        dplyr::filter(`SNP Name` %in% codes[which(codes$chr == "Y"),]$marker)

    #### SAFE STOPPING POINT
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
        # cX[gX[wh,"SNP Name"],i] <- (gX$X[wh] + gX$Y[wh])/2

        wh <- (gY[,"Sample ID"]==samples[i])
        cY[gY[wh,]$`SNP Name`,i] <- (gY$X[wh] + gY$Y[wh])/2
        # cY[gY[wh,"SNP Name"],i] <- (gY$X[wh] + gY$Y[wh])/2
    }
    if(is.null(cXint)) {
        cXint <- cX
        cYint <- cY
    } else {
        # if any columns in both, use those from second set
        cXint <- qtl2convert::cbind_smother(cXint, cX)
        cYint <- qtl2convert::cbind_smother(cYint, cY)
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

# Write all intensities
# unzip and read the data
dat <- vector("list", length(ifiles))
for(i in seq_along(ifiles)) {
  # zipfile <- tempfile()
  # download.file(ifiles[i], zipfile)
  unzipped <- unzip(ifiles[i])
  dat[[i]] <- vroom::vroom(file = unzipped, skip = 9,
               num_threads = parallel::detectCores())
  # unlink(zipfile)
}

# rbind the results together, saving selected columns
dat <- do.call("rbind", dat)[,c("SNP Name", "Sample ID", "X", "Y")]

# create matrices that are snps x samples
snps <- unique(dat[,"SNP Name"])
samples <- unique(dat[,"Sample ID"])
X <- Y <- matrix(ncol=length(samples$`Sample ID`), nrow=length(snps$`SNP Name`))
dimnames(X) <- dimnames(Y) <- list(snps$`SNP Name`, samples$`Sample ID`)
for(i in seq(along=samples$`Sample ID`)) {
  message(i, " of ", length(samples$`Sample ID`))
  tmp <- dat[dat[,"Sample ID"]==samples$`Sample ID`[i],]
  X[,samples$`Sample ID`[i]] <- tmp[,"X"]$X
  Y[,samples$`Sample ID`[i]] <- tmp[,"Y"]$Y
}

# bring together in one matrix
result <- cbind(snp=rep(snps$`SNP Name`, 2),
                channel=rep(c("X", "Y"), each=length(snps$`SNP Name`)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps$`SNP Name`), seq_along(snps$`SNP Name`)+length(snps$`SNP Name`)))),]
rownames(result) <- 1:nrow(result)

# write to fst file, maximally compressed
write_fst(result, "intensities.fst", compress=100)






