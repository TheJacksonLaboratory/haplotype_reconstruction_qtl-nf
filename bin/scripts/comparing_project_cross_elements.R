library(qtl2)
library(dplyr)
library(tidyr)

# load a successful cross
load("/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/do_oocyte/geno_probs/working_cross.RData")
oocyte_cross <- working_cross

# load the bad cross
load("/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/DO_ESC/geno_probs/working_cross.RData")
esc_cross <- working_cross
rm(working_cross)

# cross type
oocyte_cross$crosstype == esc_cross$crosstype

# check genotype encoding
chrom = c(1:length(oocyte_cross$geno))

oocyte_summary <- purrr::map2(.x = chrom,
            .y = oocyte_cross$geno, 
            .f = function(x,y){
              chr_table <- apply(y, 2, function(x) {
                length(table(x))})
              data.frame(chr_table) %>%
                dplyr::mutate(chrom = x)
              }) %>%
  Reduce(rbind,.)

oocyte_summary %>%
  dplyr::group_by(chrom,chr_table) %>%
  dplyr::count() %>%
  tidyr::pivot_wider(names_from = chr_table, values_from = n)


esc_summary <- purrr::map2(.x = chrom,
                             .y = esc_cross$geno, 
                             .f = function(x,y){
                               chr_table <- apply(y, 2, function(x) {
                                 length(table(x))})
                               data.frame(chr_table) %>%
                                 dplyr::mutate(chrom = x)
                             }) %>%
  Reduce(rbind,.)

esc_summary %>%
  dplyr::group_by(chrom,chr_table) %>%
  dplyr::count() %>%
  tidyr::pivot_wider(names_from = chr_table, values_from = n)


# genetic map comparison
oocyte_cross$gmap[[1]][names(oocyte_cross$gmap[[1]]) %in% names(esc_cross$gmap[[1]])][1:10]
esc_cross$gmap[[1]][names(esc_cross$gmap[[1]]) %in% names(oocyte_cross$gmap[[1]])][1:10]

# physical map comparison
oocyte_cross$pmap[[1]][names(oocyte_cross$pmap[[1]]) %in% names(esc_cross$pmap[[1]])][1:10]
esc_cross$pmap[[1]][names(esc_cross$pmap[[1]]) %in% names(oocyte_cross$pmap[[1]])][1:10]

# covar files
head(oocyte_cross$covar)
head(esc_cross$covar)

# cross info
oocyte_cross$cross_info
esc_cross$cross_info