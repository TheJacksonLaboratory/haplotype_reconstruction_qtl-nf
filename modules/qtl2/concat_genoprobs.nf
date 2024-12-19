process CONCAT_GENOPROBS {

  cpus 4
  memory {500.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  publishDir "${params.projectDir}/projects/${project_id}/results", pattern:"*.rds", mode:'copy'

  input:
  tuple file(crosses), val(project_id), file(excluded_samples), file(genoprobs)

  output:
  tuple file(crosses), val(project_id), file(excluded_samples), file("*_genoprobs.rds"), file("*_alleleprobs.rds"), file("*_cross.rds"), emit: concat_probs

  script:
  
  """
  echo ${genoprobs} > probs.txt
  echo ${crosses} > crosses.txt
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/concatGenoProbs.R
  mv genoprobs.rds ${project_id}_genoprobs.rds
  mv alleleprobs.rds ${project_id}_alleleprobs.rds
  mv cross.rds ${project_id}_cross.rds

  """
}
