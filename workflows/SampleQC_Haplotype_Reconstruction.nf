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
include {CONCAT_INTENSITIES} from "${projectDir}/modules/qtl2/concat_intensities"
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

    // Process FinalReport File
    GS_TO_QTL2(project_ch)
    metadata = GS_TO_QTL2.out.qtl2meta
    sampleGenos = GS_TO_QTL2.out.sampleGenos
    
    // Gather intensities by project id
    sexchr_intensities = GS_TO_QTL2.out.qtl2ints
                                    .groupTuple(by: 1)
                                    .map {it -> [it[0].flatten(), it[1]]}

    all_intensities = GS_TO_QTL2.out.qtl2intsfst.groupTuple(by: 1)
    qc_intensities = sexchr_intensities.join(all_intensities, by: [1, 1])

    // Concatenate intensities across projects for QC
    CONCAT_INTENSITIES(qc_intensities)

    // Write control file
    WRITE_CROSS(metadata, sampleGenos, consensusFiles)

    // Initial haplotype reconstruction
    GENOPROBS(WRITE_CROSS.out.cross)

    // Gather genoprobs by project id
    project_genoprobs = GENOPROBS.out.genoprobs.groupTuple(by: 1)

    // Concatenate genoprobs across projects and perform marker QC
    CONCAT_GENOPROBS(project_genoprobs)

    // Join all the hr elements
    qc_data = CONCAT_GENOPROBS.out.concat_probs.join(CONCAT_INTENSITIES.out.concat_intensities)

    // Render the QC report
    QC_REPORT(qc_data)

}
