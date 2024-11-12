process SAMPLE_MARKER_QC {

  cpus 1
  memory {50.GB * task.attempt}
  time {2.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 3

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.sample_folder}/geno_probs", pattern: "*.RData", mode:'copy'

  input:
  tuple file(cross), file(intensities)

  output:
  path("QC_1.RData"), emit: qc_data
  path("working_cross.RData"), emit: genoprobs_cross

  script:
  log.info "----- Performing Initial Sample and Marker Quality Control -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/sampleQC.R ${cross} ${intensities}

  """
}
