#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Nextflow pipeline for sample QC and haplotype reconstruction
// on genetically diverse mice

// import modules
include {COUNT_SAMPLES} from "${projectDir}/modules/general/count_samples"
include {GS_TO_QTL2} from "${projectDir}/modules/qtl2/geneseek2qtl2"
include {WRITE_CROSS} from "${projectDir}/modules/qtl2/write_cross"
include {GENOPROBS} from "${projectDir}/modules/qtl2/genoprobs"
include {CONCAT_GENOPROBS} from "${projectDir}/modules/qtl2/concat_genoprobs"
include {QC_REPORT} from "${projectDir}/modules/markdown/render_QC_markdown"

// include {SAMPLE_MARKER_QC} from "${projectDir}/modules/qtl2/sample_marker_QC"

// hold for including a help page if help if needed
// if (params.help){
//     help()
//     exit 0
//}

// hold for printing the log of selected parameters
// param_log()


// save for starting point
// check files for errors

// Make channel of consensus files (GigaMUGA)
GM_foundergenos = Channel.fromPath("${params.CCDOdataDir}/GM_foundergeno*").collect()
GM_gmaps = Channel.fromPath("${params.CCDOdataDir}/GM_gmap*").collect()
GM_pmaps = Channel.fromPath("${params.CCDOdataDir}/GM_pmap*").collect()
consensusFiles = GM_foundergenos
                    .concat(GM_gmaps)
                    .concat(GM_pmaps)
                    .flatten().collect()


// QC and Haplotype Reconstruction Workflow
workflow HR_QC {
    
    // Create channels for each FinalReport file
    project_ch = Channel.fromPath("${params.manifest}")
                        .splitCsv(header: true)
                        .map {row -> 
                        [ finalreport_file = row.finalreport_file,
                          project_id = row.project_id,
                          covar_file = row.covar_file,
                          cross_type = row.cross_type ]}
                        // *** snippet to process by project instead of by FinalReport ***
                        //.groupTuple(by: 1)
                        //.map {it -> [it[0], it[1], it[2].unique().flatten()[0], it[3].unique().flatten()[0]]}

    // Process FinalReport File
    GS_TO_QTL2(project_ch)
    intensities = GS_TO_QTL2.out.qtl2intsfst.groupTuple(by: 1)
    metadata = GS_TO_QTL2.out.qtl2meta
    sampleGenos = GS_TO_QTL2.out.sampleGenos

    // Write control file
    WRITE_CROSS(metadata, sampleGenos, consensusFiles)

    // Initial haplotype reconstruction
    GENOPROBS(WRITE_CROSS.out.cross)

    // Gather by project id
    project_genoprobs = GENOPROBS.out.genoprobs.groupTuple(by: 1)

    // Concatenate genoprobs across projects and perform marker QC
    CONCAT_GENOPROBS(project_genoprobs)





    
    
    // Perform initial sample QC
    //SAMPLE_MARKER_QC(crosses)



    // Render the QC report
    //report_data = GENOPROBS_QC.out.genoprob_qc
    //			.combine(WRITE_CROSS.out.cross)
    //			.combine(GS_TO_QTL2.out.qtl2ints)
    //			.combine(SAMPLE_MARKER_QC.out.qc_data)
    //report_data.view()
    //QC_REPORT(report_data)

}
