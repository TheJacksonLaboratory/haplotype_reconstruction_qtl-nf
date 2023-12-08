#!/usr/bin/env Rscript
library(qtl2)
library(qtl2convert)
################################################################################
# Make .json file and cross object from genotype files and GigaMUGA reference
# data.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20231208
################################################################################ 
# test_dir <- "/fastscratch/QC_HAP_outputDir/work/b2/23a342cd54730a19c9efab3e958378"
# setwd(test_dir)
args <- commandArgs(trailingOnly = TRUE)


# covar file
cat(" -Input files:\n")
covar_file <- args[1]
# covar_file <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/do_oocyte/covar_files/DO_covar_nf.csv"
print(covar_file)

# founder data
founder_ostem <- args[2]
# founder_ostem <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data"
print(founder_ostem)

# make the control file in the output directory?
sample_geno_ostem <- args[3]
# sample_geno_ostem <- '/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/do_oocyte/qtl2genos'
print(sample_geno_ostem)

# read in covariate file
covar <- read.csv(covar_file, stringsAsFactors = F)
covar$sex[covar$sex == "female"] <- "F"
covar$sex[covar$sex == "male"] <- "M"
write.csv(x = covar, file = "DO_covar_nf.csv", row.names = F)
head(covar)


# move elements to one directory for reference and for proper cross paths
system(paste0("mv GM*.csv ",sample_geno_ostem))
system(paste0("mv geno*.csv ",sample_geno_ostem))
system(paste0("mv DO_covar_nf.csv ",sample_geno_ostem))

# Write control file
chr <- c(1:19, "X")
cat(" -Writing control file\n")
# for now stick with "do" cross type, eventually take a crosstype param
# have to paste getwd() in front of the control file output name so that
# the files in the channel are recognized properly by read_cross2

if(unique(covar$sex) == "female"){
  cat(" Note: only females present in cross.")
  qtl2::write_control_file(output_file = file.path(sample_geno_ostem,"QC_HAP.json"),
                           crosstype="do",
                           # description="QC_HAP",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="DO_covar_nf.csv",
                           crossinfo_covar=colnames(covar)[!colnames(covar) %in% "id"],
                           sex_covar="sex",
                           sex_codes=c(F="female"),
                           overwrite = T)
} else if(unique(covar$sex) == "male"){
  cat(" Note: only males present in cross.")
  qtl2::write_control_file(output_file = file.path(sample_geno_ostem,"QC_HAP.json"),
                           crosstype="do",
                           # description="QC_HAP",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="DO_covar_nf.csv",
                           crossinfo_covar=colnames(covar)[!colnames(covar) %in% "id"],
                           sex_covar="sex",
                           sex_codes=c(M="male"),
                           overwrite = T)
} else {
  qtl2::write_control_file(output_file = file.path(sample_geno_ostem,"QC_HAP.json"),
                           crosstype="do",
                           # description="QC_HAP",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="DO_covar_nf.csv",
                           crossinfo_covar=colnames(covar)[!colnames(covar) %in% "id"],
                           sex_covar="sex",
                           sex_codes=c(F="female", M="male"),
                           overwrite = T)
}


cat(" -Reading control file\n")
# read cross
cross <- qtl2::read_cross2(file = file.path(sample_geno_ostem,"QC_HAP.json"), quiet = F)
save(cross, file = "preQC_cross.RData")

