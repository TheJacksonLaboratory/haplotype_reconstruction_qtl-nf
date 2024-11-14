process WRITE_CROSS {

  cpus 16
  memory 50.GB
  time '01:00:00'
  errorStrategy { task.exitStatus == 138..143 ? 'retry' : 'terminate' }
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${projectDir}/results/${project_id}/qtl2genos", pattern: "*.RData", mode:'copy'

  input:
  tuple file(covar), val(project_id), path(covar_file), val(cross_type)
  path(sampleGenos)
  path(consensusFiles)

  output:
  tuple path("*.RData"), val(project_id), emit: cross

  script:
  log.info "----- Making Control Files, R/qtl2 Cross Object: Project ${project_id} -----"

  """
  echo ${consensusFiles}

  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/writeControlFile.R ${projectDir}/results/${project_id}/qtl2genos
  """
}
