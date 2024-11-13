#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Nextflow pipeline for sample QC and haplotype reconstruction
// on genetically diverse mice

// import modules
include {GS_TO_QTL2} from "${projectDir}/modules/qtl2/geneseek2qtl2"
include {WRITE_CROSS} from "${projectDir}/modules/qtl2/write_cross"
include {SAMPLE_MARKER_QC} from "${projectDir}/modules/qtl2/sample_marker_QC"
include {GENOPROBS_QC} from "${projectDir}/modules/qtl2/genoprobs_qc"
include {QC_REPORT} from "${projectDir}/modules/markdown/render_QC_markdown"

// hold for including a help page if help if needed
// if (params.help){
//     help()
//     exit 0
//}

// hold for printing the log of selected parameters
// param_log()


// save for starting point
// check files for errors



//FinalReports = Channel.fromPath("${params.sample_folder}/neogen_finalreports/*FinalReport*").collect()
GM_foundergenos = Channel.fromPath("${params.CCDOdataDir}/GM_foundergeno*").collect()
GM_gmaps = Channel.fromPath("${params.CCDOdataDir}/GM_gmap*").collect()
GM_pmaps = Channel.fromPath("${params.CCDOdataDir}/GM_pmap*").collect()

// QC and Haplotype Reconstruction Workflow
workflow QC_HAP {
    

    // Create channels
    project_ch = Channel.fromPath("${params.manifest}")
                    .splitCsv(header: true)
                    .map {row -> 
                            [ finalreport_file = row.finalreport_file.toString(),
                            project_id = row.project_id.toString(),
                            covar_file = row.covar_file.toString(),
                            cross_type = row.cross_type.toString() ]}
                    .groupTuple(by: 1)
                    .map {it -> [it[0], it[1], it[2].unique().flatten()[0], it[3].unique().flatten()[0]]}
    project_ch.view()

    // Process FinalReport File
    GS_TO_QTL2(FinalReports)

    // Write control file
    //cross_elements = GS_TO_QTL2.out.qtl2genos
    //					.mix(GM_foundergenos,GM_gmaps,GM_pmaps)
    //					.flatten()
    //					.collect()

    //cross_elements.view()
    //WRITE_CROSS(cross_elements)


    // Perform initial sample QC
    //sample_QC_files = WRITE_CROSS.out.cross
    //					.combine(GS_TO_QTL2.out.qtl2intsfst)

    //SAMPLE_MARKER_QC(sample_QC_files)

    // Initial haplotype reconstruction for genotyping errors and crossover estimation
    //GENOPROBS_QC(SAMPLE_MARKER_QC.out.genoprobs_cross)


    // Render the QC report
    //report_data = GENOPROBS_QC.out.genoprob_qc
    //			.combine(WRITE_CROSS.out.cross)
    //			.combine(GS_TO_QTL2.out.qtl2ints)
    //			.combine(SAMPLE_MARKER_QC.out.qc_data)
    //report_data.view()
    //QC_REPORT(report_data)

}
