#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=haplotype_reconstruction-nf
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow main.nf \
        --workflow SampleQC_Haplotype_Reconstruction  \
        --sample_folder 'projects/DO_ESC' \
        --pubdir '/flashscratch/widmas/QC_HAP_outputDir' \
        -w '/fastscratch/QC_HAP_outputDir/work' \
        --crossType 'do' \
        --covar 'projects/DO_ESC/covar_files/DO_ESC_covar.csv' \
        --comment "This script will perform sample QC and haplotype reconstruction on genetically diverse mouse samples" \
        -resume
