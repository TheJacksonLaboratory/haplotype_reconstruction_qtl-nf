process GENOPROBS {

  cpus 4
  memory 360.GB
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple path(cross), val(project_id), file(excluded_samples)

  output:
  tuple file(cross), val(project_id), file(excluded_samples), file("*.rds"), emit: genoprobs

  script:
  
  """
  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  echo \$hash
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/calcGenoProbs.R
  mv pr_36state.rds pr_36state_\$hash.rds

  """
}
