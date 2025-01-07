process GS_TO_QTL2 {
  
  cpus 2
  time 1.hour
  memory 50.GB
  errorStrategy {(task.exitStatus == 1) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}.\n Please check ${projectDir}/results/${project_id}/logs.\n Also, please verify that sample names in metadata match those expected in FinalReport file(s).\n\n"; return 'ignore'}.call() : 'finish'}

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple path(finalreport_files), val(project_id), path(covar_file), val(cross_type)

  output:
  tuple path("*geno*.csv"), file("excluded_samples.csv"), emit: sampleGenos
  tuple file("covar.csv"), val(project_id), val(cross_type), emit: qtl2meta
  tuple path("*int.csv"), val(project_id), emit: qtl2ints
  tuple path("*.fst"), val(project_id), emit: qtl2intsfst


  script:
  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/geneseek2qtl2_nf.R \
	${params.CCDOdataDir} \
	${covar_file} \
	${finalreport_files} \
  ${params.max_pct_missing}

  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  mv intensities.fst intensities_\${hash}.fst
  mv chrYint.csv chrY_\${hash}_int.csv
  mv chrXint.csv chrX_\${hash}_int.csv
  """
}
