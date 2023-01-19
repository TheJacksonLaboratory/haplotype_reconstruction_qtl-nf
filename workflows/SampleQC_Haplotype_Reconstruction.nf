#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for sample QC and haplotype reconstruction
// on genetically diverse mice
// General flow:
// 1) 

// import modules
// include {help} from "${projectDir}/bin/help/wgs.nf"
// include {param_log} from "${projectDir}/bin/log/stitch.nf"
include {GS_TO_QTL2} from "${projectDir}/modules/qtl2/geneseek2qtl2"
include {WRITE_CONTROL_FILE} from "${projectDir}/modules/qtl2/writeControlFile"
include {SAMPLE_QC} from "${projectDir}/modules/qtl2/sampleQC"
include {CALC_GENO_PROBS} from "${projectDir}/modules/qtl2/calcGenoProbs"
include {PLOT_GENO_PROBS} from "${projectDir}/modules/qtl2/plotGenoProbs"
include {RENDER_MARKDOWN} from "${projectDir}/modules/qtl2/plotGenoProbs"



// hold for including a help page if help if needed
// if (params.help){
//     help()
//     exit 0
//}

// hold for printing the log of selected parameters
// param_log()


// save for starting point
// if any channel is empty give error message and exit
// read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

chrs = Channel.of(1..19,'X')

// QC and Haplotype Reconstruction Workflow
workflow QC_HAP {}
