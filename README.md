# haplotype_reconstruction_qtl-nf: A slim Nextflow pipeline for mouse cross haplotype reconstruction and quality control

JAX users are required to have access to the Sumner cluster, and to have Nextflow installed in their home directory. Any setup for external users will require additional support, and those wishing to share these workflows are encouraged to contact the maintainers of this repository.

This pipeline is implemented using [Nextflow](https://www.nextflow.io/), a scalable, reproducible, and increasingly common language used in the development and maintenance of bioinformatics workflows. The modular nature of the workflow is enabled by software containers, such as [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity), with all the software requirements for executing each step. Specific combinations and versions of software are specified in each container making analyses perfectly reproducible over time as long as the source data is unchanged.

## Execution:

Clone the repository using the standard procedure. On the JAX HPC, from within the cloned `haplotype_reconstruction_qtl-nf` directory:

``` bash
sbatch run_scripts/example_run_script.sh
```

## Overview:

The pipeline reads in the raw genotypes from GigaMUGA FinalReport files and makes them into files amenable to analysis using R/qtl2. These include cross files, genotype probabilities, allele probabilities, imputed genotypes from probabilities (`maxmarg` output).

Files used for sample and genotype quality control are also generated, such as inferred genotyping errors, poorly performing markers, and a markdown document outlining results from sex checks and calculations of sample duplication.

```mermaid
flowchart TD
    p0((FinalReport.zip))
    p1((R/qtl2 Covar File))
    p2((GigaMUGA Reference Files))
    p3((Phenotype Table))
    p4[GS_TO_QTL2]:::process
    p5[WRITE_CROSS]:::process
    p6[GENOPROBS]:::process
    p7[CONCAT_GENOPROBS]:::process
    p8[CONCAT_INTENSITIES]:::process
    p9[PHENOTYPE_QC]:::process
    p10[PERFORM_PERMUTATIONS]:::process
    p11[MAP_QTL]:::process
    p12[PROCESS_QTL]:::process
    o1((Sex Chromosome Marker Intensities)):::output
    o2((All Marker Intensities)):::output
    o3((Excluded File List)):::output
    o4((36-state/Genotype Probabilities)):::output
    o5((8-state/Allele Probabilities)):::output
    o6((Cross Object)):::output
    o7((Imputed Genotype States)):::output
    o8((Genotyping Error LOD Scores)):::output
    o9((Sample Quality Control Flag Summary)):::output
    o10((Bad Markers List)):::output
    o11((Sample Quality Control Summary Markdown)):::output
    o12((Phenotype Quality Control Summary Markdown)):::output
    o13((Processed Phenotypes)):::output
    o14((Kinship Matrix)):::output
    o15((Permutation Threshold Summary)):::output
    o16((QTL Table)):::output
    o17((scan1 Output)):::output
    o18((scan1blup Output)):::output
    o19((QTL LOD Plots)):::output
    o20((QTL Effect Plots)):::output
    o21((Cross Object)):::output

    p0 --> p4
    p1 --> p4

    p4 --> p5
    p4 --> p8
    p4 --> o1
    p4 --> o2

    p2 --> p4
    p4 --> p6
    p5 --> o21
    p6 --> p7

    p7 --> o3
    p7 --> o4
    p7 --> o5
    p7 --> o6
    p7 --> o7
    p7 --> o8

    p8 --> o9
    p8 --> o10
    p8 --> o11

    p3 --> p9
    p9 --> o12
    p9 --> o13
    p9 --> p10

    p5 --> p10
    p7 --> p10
    p10 --> o14
    p10 --> o15
    p10 --> p11

    p11 --> o16
    p11 --> o17
    p11 --> o19
    p11 --> p12

    p12 --> o18
    p12 --> o20

classDef output fill:#99e4ff,stroke:#000000,stroke-width:5px,color:#000000
classDef process fill:#00A2DC,stroke:#000000,stroke-width:2px,color:#000000
```

The run script `example_run_script.sh` specifies only one user--generated comma-separated sample manifest with four named columns: **finalreport_file**, **project_id**, **covar_file**, and **cross_type** (see README within `sample_sheets` subdirectory).
