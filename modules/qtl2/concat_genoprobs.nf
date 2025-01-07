process CONCAT_GENOPROBS {

  cpus 6
  memory {360.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_genoprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_alleleprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_cross.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_maxmarg.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_genotyping_errors.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_excluded_samples.csv", mode:'copy'

  input:
  tuple file(crosses), val(project_id), file(excluded_samples), file(genoprobs)

  output:
  tuple val(project_id), file("*_excluded_samples.csv"), file("*_genoprobs.rds"), file("*_alleleprobs.rds"), file("*_cross.rds"), file("*_maxmarg.rds"), file("*_genotyping_errors.rds"), emit: concat_probs

  script:
  
  """
  echo ${genoprobs} > probs.txt
  echo ${crosses} > crosses.txt
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/concatGenoProbs.R
  mv genoprobs.rds ${project_id}_genoprobs.rds
  mv alleleprobs.rds ${project_id}_alleleprobs.rds
  mv cross.rds ${project_id}_cross.rds
  mv maxmarg.rds ${project_id}_maxmarg.rds
  mv genotyping_errors.rds ${project_id}_genotyping_errors.rds
  mv excluded_samples.csv ${project_id}_excluded_samples.csv

  """
}
