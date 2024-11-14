#!/usr/bin/env Rscript
# load required packages
library(data.table)
library(qtl2convert)
library(vroom)
library(parallel)
library(dplyr)
library(magrittr)
library(fst)
args <- commandArgs(trailingOnly = TRUE)
cat(paste0(" -- R/qtl2 version: ",qtl2::qtl2version(),"\n"))
cat(paste0(" -- R/qtl2convert version: ",packageVersion("qtl2convert"),"\n"))

# GigaMUGA reference data directory
GM_ref_dir <- args[1]

# metadata file indicating which mice should be retained from array files
metadata_path <- args[2]

# FinalReport files
ifiles <- args[3]

# testing:
# GM_ref_dir <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data"
# metadata_path <- "/flashscratch/widmas/QC_HAP_outputDir/work/56/68e36632c9cfe104234dbe6b2db36c/DO_ESC_covar.csv"
# ifiles <- "/flashscratch/widmas/QC_HAP_outputDir/work/56/68e36632c9cfe104234dbe6b2db36c/finalreportlist.txt"

# allele codes
codefile <- file.path(GM_ref_dir,"GM_allelecodes.csv")

cat("Reading covar file")
metadata <- read.csv(metadata_path)

# input files with GigaMUGA genotypes
# read in the FinalReport file table from nextflow
finalreports <- read.table(ifiles)
ifiles <- c(as.matrix(t(finalreports)))

# read genotype codes
codes <- data.table::fread(codefile, skip = 3, data.table = F)

full_geno <- NULL
cXint <- cYint <- NULL
excluded_samples <- vector("list", length(ifiles))
for(ifile in ifiles) {
  cat(" -File:", ifile, "\n")
  rezip <- FALSE
  if(!file.exists(ifile)) {
    cat(" -Unzipping file\n")
    system(paste("gunzip", ifile))
    rezip <- TRUE
  }
  
  cat(" -Reading data\n")
  g <- vroom::vroom(file = ifile, skip = 9,
                    num_threads = parallel::detectCores()/1.2)
  g %<>%
    dplyr::filter(`SNP Name` %in% codes$marker,
                  `Sample ID` %in% metadata$id)
  
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
  
  # find samples with lots of missing genos off the bat
  missing_geno_samples <- which(colSums(geno == "-")/nrow(geno) > 0.10)
  cat(paste0(length(missing_geno_samples)," samples missing more than 10% of genotypes \n"))
  if(length(missing_geno_samples) > 0){
    geno <- geno[,-missing_geno_samples]
    excluded_samples[[which(ifile == ifiles)]] <- missing_geno_samples
  }
  
  
  if(is.null(full_geno)) {
    full_geno <- geno
  } else {
    # if any columns in both, use those from second set
    full_geno <- cbind(full_geno, geno)
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
    cXint <- cbind(cXint, cX)
    cYint <- cbind(cYint, cY)
  }
  
  if(rezip) {
    cat(" -Rezipping file\n")
    system(paste("gzip", ifile))
  }
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(df = cbind(marker=rownames(cXint), cXint), 
                       filename = "chrXint.csv", 
                       comment = "X chr intensities",
                       overwrite=TRUE)
qtl2convert::write2csv(df = cbind(marker=rownames(cYint), cYint),
                       filename = "chrYint.csv",
                       comment ="Y chr intensities",
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chrom in c(1:19,"X","Y","M")) {
  mar <- codes %>%
    dplyr::filter(chr == chrom) %>%
    dplyr::select(marker)
  g <- full_geno[rownames(full_geno) %in% mar$marker,]
  qtl2convert::write2csv(df = cbind(marker=rownames(g), g),
                         filename = paste0("geno", chrom, ".csv"),
                         comment = paste0("genotypes for chr ", chrom),
                         overwrite=TRUE)
}


# Write all intensities
# unzip and read the data
dat <- vector("list", length(ifiles))
for(i in seq_along(ifiles)) {
  unzipped <- unzip(ifiles[i])
  pre_dat <- vroom::vroom(file = unzipped, skip = 9,
                          num_threads = parallel::detectCores()/1.2)
  ex_samples <- excluded_samples[[i]]
  if(length(ex_samples) > 0){
    dat[[i]] <- pre_dat %>%
      dplyr::filter(`SNP Name` %in% rownames(full_geno),
                    `Sample ID` %in% colnames(full_geno),
                    !`Sample ID` %in% names(ex_samples))
  } else {
    dat[[i]] <- pre_dat %>%
      dplyr::filter(`SNP Name` %in% rownames(full_geno),
                    `Sample ID` %in% colnames(full_geno))
  }
}

# rbind the results together, saving selected columns
dat <- do.call("rbind", dat)[,c("SNP Name", "Sample ID", "X", "Y")]

# create matrices that are snps x samples
snps <- unique(dat$`SNP Name`)
samples <- unique(dat$`Sample ID`)
X <- Y <- matrix(ncol=length(samples), 
                 nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
  cat(i, " of ", length(samples))
  tmp <- dat[which(dat$`Sample ID`==samples[i]),]
  X[,samples[i]] <- tmp$X
  Y[,samples[i]] <- tmp$Y
}

# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), 
                            each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps),
                                    seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

# write to fst file, maximally compressed
write_fst(result, "intensities.fst", compress=100)

# write the covariate file with a generic name
write.csv(metadata, "covar.csv", quote = F, row.names = F)

