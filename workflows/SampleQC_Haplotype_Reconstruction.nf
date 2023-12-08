#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for sample QC and haplotype reconstruction
// on genetically diverse mice

// import modules
// include {help} from "${projectDir}/bin/help/wgs.nf"
// include {param_log} from "${projectDir}/bin/log/stitch.nf"

include {GS_TO_QTL2} from "${projectDir}/modules/qtl2/geneseek2qtl2"
include {WRITE_CROSS} from "${projectDir}/modules/qtl2/write_cross"
include {SAMPLE_MARKER_QC} from "${projectDir}/modules/qtl2/sample_marker_QC"

// hold for including a help page if help if needed
// if (params.help){
//     help()
//     exit 0
//}

// hold for printing the log of selected parameters
// param_log()


// save for starting point
// check files for errors

chrs = Channel.of(1..19,'X')
FinalReports = Channel.fromPath("${params.sample_folder}/neogen_finalreports/*FinalReport*").collect()
GM_foundergenos = Channel.fromPath("${params.CCDOdataDir}/GM_foundergeno*").collect()
GM_gmaps = Channel.fromPath("${params.CCDOdataDir}/GM_gmap*").collect()
GM_pmaps = Channel.fromPath("${params.CCDOdataDir}/GM_pmap*").collect()

// QC and Haplotype Reconstruction Workflow
workflow QC_HAP {
    
    // Process FinalReport File
    GS_TO_QTL2(FinalReports)

    // Write control file
    cross_elements = GS_TO_QTL2.out.qtl2genos
					.mix(GM_foundergenos,GM_gmaps,GM_pmaps)
					.flatten()
					.collect()
    WRITE_CROSS(cross_elements)


    // Perform initial sample QC
    sample_QC_files = WRITE_CROSS.out.cross.combine(GS_TO_QTL2.out.qtl2intsfst)
    SAMPLE_MARKER_QC(sample_QC_files)

}
