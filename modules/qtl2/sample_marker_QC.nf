process SAMPLE_MARKER_QC {

  cpus 1
  memory {50.GB * task.attempt}
  time {1.hour * task.attempt}
  errorStrategy { task.exitStatus == 138..143 ? 'retry' : 'terminate' }
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${projectDir}/results/${project_id}/QC", pattern: "*.RData", mode:'copy'

  input:
  tuple val(project_id), path(preQCdata), val(nsamples), file(intensities)

  output:
  tuple path("QC_1.RData"), val(project_id), emit: qc_data
  tuple path("working_cross.RData"), val(project_id), emit: genoprobs_cross

  script:
  log.info "----- Performing Initial Sample and Marker Quality Control: Project ${project_id} -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/sampleQC.R ${preQCdata} ${intensities}

  """
}
