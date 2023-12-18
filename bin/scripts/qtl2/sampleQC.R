#!/usr/bin/env Rscript
library(qtl2)
library(dplyr)
library(ggplot2)
library(fst)
################################################################################
# Perform marker and sample QC using cross object and intensities upstream in 
# haplotype reconstruction pipeline.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20231208
################################################################################ 
# test_dir <- "/fastscratch/QC_HAP_outputDir/work/b8/fdfd8ce7bdd497445041b0ac27a8b7"
# setwd(test_dir)
args <- commandArgs(trailingOnly = TRUE)

# import cross object
cross <- args[1]
# cross <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/do_oocyte/geno_probs/cross.RData"
load(cross)

# import intensities
intensities <- args[2]
# intensities <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/projects/do_oocyte/qtl2genos/intensities.fst"
print(intensities)


# Reordering genotypes so that most common allele in founders is first
for(chr in seq_along(cross$founder_geno)) {
  fg <- cross$founder_geno[[chr]]
  g <- cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  cross$founder_geno[[chr]] <- fg
  cross$geno[[chr]] <- g
}

# percent missing genotypes - DO
percent_missing <- qtl2::n_missing(cross, "ind", "prop")*100

# Sample Duplicates
cg <- compare_geno(cross, cores=0)

# Sex checks
# Reading in all probe intensities
int <- fst::read.fst(intensities)
int <- int[seq(1, nrow(int), by=2),-(1:2)] + int[-seq(1, nrow(int), by=2),-(1:2)]
int <- int[,which(colnames(int) %in% qtl2::ind_ids(cross))]

## Genotype Frequencies
g <- do.call("cbind", cross$geno[1:19])
fg <- do.call("cbind", cross$founder_geno[1:19])

# find markers with missing genotypes in samples
g <- g[,colSums(fg==0)==0]

# find markers with missing genotypes in the founders
fg <- fg[,colSums(fg==0)==0]
fgn <- colSums(fg==3)

# save objects to send to markdown
save(fg, # markers with missing genotypes in the founders
     g,  # markers with missing genotypes in a sample
     cg, # sample duplicates
     int, # marker intensities
     percent_missing, # percent of missing genotypes per sample
     file = "QC_1.RData")

# gf_ind <- vector("list", 4)
# for(i in 1:4) {
#   gf_ind[[i]] <- t(apply(g[,fgn==i], 1, function(a) table(factor(a, 1:3))/sum(a != 0)))
# }

# Interactive plot of DO array intensities per sample
# n <- names(sort(percent_missing, decreasing=TRUE))
# iboxplot <- iboxplot(log10(t(int)+1),
#                      orderByMedian=TRUE,
#                      chartOpts=list(ylab="log10(SNP intensity + 1)"))

# # plot genotype frequencies
# par(mfrow=c(2,2), mar=c(0.6, 0.6, 2.6, 0.6))
# for(i in 1:4) {
#   triplot(c("AA", "AB", "BB"), main=paste0("MAF = ", i, "/8"))
#   tripoints(gf_ind[[i]], pch=21, bg="lightblue")
#   tripoints(c((1-i/8)^2, 2*i/8*(1-i/8), (i/8)^2), pch=21, bg="violetred")
# 
#   if(i>=3) { # label mouse with lowest het
#     wh <- which(gf_ind[[i]][,2] == min(gf_ind[[i]][,2]))
#     tritext(gf_ind[[i]][wh,,drop=FALSE] + c(0.02, -0.02, 0),
#             names(wh), adj=c(0, 1))
#   }
# 
#   # label other mice
#   if(i==1) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.3]
#   }
#   else if(i==2) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.48]
#   }
#   else if(i==3) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.51]
#   }
#   else if(i==4) {
#     lab <- rownames(gf_ind[[i]])[gf_ind[[i]][,2]>0.6]
#   }
# 
#   for(ind in lab) {
#     if(grepl("^F", ind) && i != 3) {
#       tritext(gf_ind[[i]][ind,,drop=FALSE] + c(-0.01, 0, +0.01), ind, adj=c(1,0.5))
#     } else {
#       tritext(gf_ind[[i]][ind,,drop=FALSE] + c(0.01, 0, -0.01), ind, adj=c(0,0.5))
#     }
#   }
# }

# xint <- xint[,colnames(xint)[which(colnames(xint) %in% metadata$sample)]]
# yint <- yint[,colnames(yint)[which(colnames(yint) %in% metadata$sample)]]
# sex <- cross$is_female
# sex <- dplyr::if_else(condition = sex == TRUE, true = "F", false = "M")
# names(sex) <- names(cross$is_female)

# if(length(levels(as.factor(sex))) == 1){
#   print("Skipping t-test; only 1 level to t-test")
#   DOxint_ave <- colMeans(DOxint, na.rm=TRUE)
#   DOyint_ave <- colMeans(DOyint, na.rm=TRUE)
#   xyints <- data.frame(DOxint_ave, DOyint_ave)
#   
#   # test
#   xyints <- cbind(xyints, sex)
#   labels <- paste0(rownames(xyints), " (", round(percent_missing), "%)")
#   DOsexcheck_plot <- ggplot(data = xyints, mapping = aes(x = DOxint_ave, 
#                                                          y = DOyint_ave, 
#                                                          fill = sex, 
#                                                          label = labels)) +
#     theme_bw() + 
#     geom_point(shape = 21, size = 4) + 
#     scale_fill_manual(values = c("green","purple")) + 
#     labs(x = "Average X chr intensity",
#          y = "Average Y chr intensity")
#   
#   plotly::ggplotly(DOsexcheck_plot, tooltip = "label")
# } else {
#   
#   print("Testing for uninformative markers and filtering those out")
#   x_pval <- apply(DOxint, 1, function(a) t.test(a ~ sex)$p.value)
#   y_pval <- apply(DOyint, 1, function(a) t.test(a ~ sex)$p.value)
#   DOxint_ave <- colMeans(DOxint[x_pval < 0.05/length(x_pval),], na.rm=TRUE)
#   DOyint_ave <- colMeans(DOyint[y_pval < 0.05/length(y_pval),], na.rm=TRUE)
#   
#   ## Plotting sex chromosome intensities to verify sexes
#   xyints <- data.frame(DOxint_ave, DOyint_ave) %>%
#     dplyr::mutate(sample = rownames(.))
#   rownames(xyints) <- NULL
#   labels <- paste0(names(DOxint_ave), " (", round(percent_missing), "%)")
#   DOsexcheck_plot <- ggplot(data = xyints, mapping = aes(x = DOxint_ave, 
#                                                          y = DOyint_ave, 
#                                                          fill = Sex, 
#                                                          label = labels)) +
#     theme_bw() + 
#     geom_point(shape = 21, size = 4) + 
#     scale_fill_manual(values = c("green","purple")) + 
#     labs(x = "Average X chr intensity",
#          y = "Average Y chr intensity")
#   plotly::ggplotly(DOsexcheck_plot, tooltip = "label")
# }