#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=hr-nf
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 36:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.out

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow main.nf \
        --workflow SampleQC_Haplotype_Reconstruction \
        --manifest 'sample_sheets/20250102_hr-nf_manifest.csv' \
        -w '/flashscratch/widmas/HR_QC_outputDir/work' \
        --comment "This script will perform sample QC and haplotype reconstruction on genetically diverse mouse samples" \
        -with-dag HR_QC_flow.html \
        -resume

