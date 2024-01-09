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
test_dir <- "/fastscratch/QC_HAP_outputDir/work/4a/a579758956a5926bcb90d4a2908619"

# file containing allele codes for GigaMUGA data
#   - from GM_processed_files.zip, https://doi.org/10.6084/m9.figshare.5404759
# codefile <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data/GM_allelecodes.csv"
codefile <- args[1]

# input files with GigaMUGA genotypes
#  - can be a single file or a vector of multiple files
#  - if samples appear in multiple files, the genotypes in later files
#    will be used in place of genotypes in earlier files
#  - files can be gzipped (".gz" extension)
# ifiles <- file.path(test_dir,"finalreportlist.txt")
ifiles <- args[2]

# metadata
metadata <- args[3]
# metadata <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/covar_files/DO_ESC_covar.csv"
meta <- read.csv(metadata)

# read genotype codes
codes <- data.table::fread(codefile, skip = 3)

full_geno <- NULL
cXint <- cYint <- NULL

finalreports <- read.table(ifiles)
ifiles <- c(as.matrix(t(finalreports)))
# ifiles <- paste(test_dir,c(as.matrix(t(finalreports))), sep = "/")

for(ifile in ifiles) {
    cat(" -File:", ifile, "\n")
    # rezip <- FALSE
    # if(!file.exists(ifile)) {
    #     cat(" -Unzipping file\n")
    #     system(paste("gunzip", ifile))
    #     rezip <- TRUE
    # }

    cat(" -Reading data\n")
    
    # benchmarking reading performance
    # microbenchmark::microbenchmark(unit = "milliseconds",
    #   g <- data.table::fread(paste0('unzip -p ',ifile),
    #                          skip = 9, 
    #                          nrows = 100),
    #   g <- vroom::vroom(file = ifile, skip = 9,
    #                     n_max = 100,
    #                     num_threads = parallel::detectCores()/2)
    # )
    # g <- data.table::fread(paste0('unzip -p ',ifile),
    #                        skip = 9, 
    #                        nrows = 100)
    g <- vroom::vroom(file = ifile, skip = 9,
                      # n_max = 100000,
                      num_threads = parallel::detectCores()/2)
    
    
    meta_sample_names <- unique(meta$id)
    sample_names <- unique(g$`Sample ID`)
    
    name_check <- Reduce(rbind, 
           lapply(meta_sample_names, 
                  function(x) data.frame(x, length(grep(x = sample_names, 
                                                        pattern = x)))))
    colnames(name_check) <- c("sample","matches")
    
    # test if *some* samples from metadata are in the finalreport file
    
    if(!1 %in% unique(name_check$matches)){
      next
    }
    
    if(0 %in% unique(name_check$matches)){
      cat(paste0(" NOTE: Some samples in supplied metadata not present in", ifile,"\n"))
      cat(paste0(" Filtering", ifile," to matching samples\n"))
      
      sample_matches <- name_check[which(name_check$matches != 0),1]
      
      
      g %<>%
        dplyr::filter(`SNP Name` %in% codes$marker,
                      `Sample ID` %in% name_check$sample)
      cat(paste0("Matching samples:\n"))
      print(unique(g$`Sample ID`))
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
        # if any columns in both, use those from second set
        full_geno <- qtl2convert::cbind_smother(full_geno, geno)
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

    # if(rezip) {
    #     cat(" -Rezipping file\n")
    #     system(paste("gzip", ifile))
    # }
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






