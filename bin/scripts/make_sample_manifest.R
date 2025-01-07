#!/usr/bin/env Rscript

################################################################################
# Make a .csv manifest for haplotype reconstruction from MUGA-series files.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241112
################################################################################

library(tidyverse)

# hr-nf directory
hr_nf_dir <- '/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf'

# output directory
out_dir <- file.path(hr_nf_dir,"sample_sheets")

# today
today <- gsub("-","",Sys.Date())

# neogen files
finalreport_file <- list.files(hr_nf_dir, pattern = "FinalReport", recursive = T, full.names = T)

# project ids
project_id <- unlist(lapply(finalreport_file, function(x) strsplit(x,"/")[[1]][[9]]))

# read in CSNA metadata
CSNA_metadata <- read.csv(file.path(hr_nf_dir, "projects/CSNA/CSNA_metadata.csv"))
trimmed_CSNA_metadata <- CSNA_metadata %>%
  dplyr::select(Sample.ID, Sex, DO.Generation, chrM) %>%
  dplyr::rename(id = Sample.ID,
                sex = Sex,
                gen = DO.Generation)
write.csv(trimmed_CSNA_metadata, file.path(hr_nf_dir, "projects/CSNA/csna_covar.csv"), quote = F, row.names = F)

# read in Attie metadata
attie_metadata_files <- list.files(file.path(hr_nf_dir,"projects/Attie"), pattern = "inventory", full.names = T)
attie_metadata <- lapply(attie_metadata_files, function(x) suppressMessages(readr::read_csv(x, trim_ws = T))) %>%
  Reduce(dplyr::bind_rows,.)
trimmed_attie_metadata <- attie_metadata %>%
  dplyr::select(`Unique Sample ID`,Sex,`DO Generation`,`Animal ID`) %>%
  dplyr::rename(id = `Unique Sample ID`,
                sex = Sex,
                gen = `DO Generation`)
write.csv(trimmed_attie_metadata, file.path(hr_nf_dir, "projects/Attie/covar_files/attie_covar.csv"), quote = F, row.names = F)

# read in Beamer metadata
beamer_metadata <- read.csv(file.path(hr_nf_dir, "projects/Beamer/tufts_animal_inventory.csv"), tryLogical = F)
trimmed_beamer_metadata <- beamer_metadata %>%
  dplyr::select(Animal.ID, Sex, DO.Generation, Strain) %>%
  dplyr::rename(id = Animal.ID,
                sex = Sex,
                gen = DO.Generation)
write.csv(trimmed_beamer_metadata, file.path(hr_nf_dir, "projects/Beamer/covar_files/beamer_covar.csv"), quote = F, row.names = F)

# read in Baker metadata
baker_metadata <- read.csv("/projects/compsci/vmp/USERS/widmas/Baker_DO_scRNAseq/DO_sampleIDs_update.csv", tryLogical = F)
wrangled_baker_metadata  <- baker_metadata %>%
  dplyr::rename(id = gigamuga_sampleID,
                gen = generation) %>%
  dplyr::select(id, sex, gen, everything())
write.csv(wrangled_baker_metadata, file.path(hr_nf_dir, "projects/Baker/covar_files/baker_covar.csv"), quote = F, row.names = F)
# Make total manifest

# covar files
covar_file <- unlist(lapply(project_id, function(x){
  list.files(file.path(hr_nf_dir, "projects", x, "covar_files"), pattern = "covar", recursive = T, full.names = T)
}))

# cross type
cross_type <- rep("do", length(covar_file))

# make .csv manifest
manifest <- data.frame(finalreport_file, project_id, covar_file, cross_type)

# write manifest
write.csv(manifest, file = file.path(out_dir,paste0(today,"_hr-nf_manifest.csv")), quote = F, row.names = F)

