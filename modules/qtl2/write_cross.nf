process WRITE_CROSS {

  memory 50.GB
  time 1.hour
  errorStrategy 'retry'
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple file(covar), val(project_id), val(cross_type)
  tuple path(sampleGenos), file(excluded_samples)
  path(consensusFiles)

  output:
  tuple path("*.rds"), val(project_id), file("excluded_samples_*"), emit: cross

  script:
  """
  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  echo \$hash
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/writeControlFile.R
  mv preQC_cross.rds preQC_cross_\${hash}.rds
  mv excluded_samples.csv excluded_samples_\${hash}.csv
  """
}
