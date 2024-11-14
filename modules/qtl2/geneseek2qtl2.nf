process GS_TO_QTL2 {

  cpus 16
  memory {50.GB * task.attempt}
  time {1.hour * task.attempt}
  errorStrategy { task.exitStatus == 138..143 ? 'retry' : 'terminate' }
  maxRetries 1
  
  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${projectDir}/results/${project_id}/intensities", pattern: "*int.csv", mode:'copy'
  publishDir "${projectDir}/results/${project_id}/intensities", pattern: "*.fst", mode:'copy'

  input:
  tuple path(finalreport_files), val(project_id), path(covar_file), val(cross_type)

  output:
  path("*geno*.csv"), emit: sampleGenos
  tuple file("covar.csv"), val(project_id), path(covar_file), val(cross_type), emit: qtl2meta
  tuple path("*int.csv"), val(project_id), path(covar_file), val(cross_type), emit: qtl2ints
  tuple path("*.fst"), val(project_id), path(covar_file), val(cross_type), emit: qtl2intsfst
  

  script:
  log.info "----- Processing FinalReport Files: Project ${project_id} -----"

  """
  echo ${finalreport_files} > finalreportlist.txt

  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/geneseek2qtl2_nf.R \
	${params.CCDOdataDir} \
	${covar_file} \
	finalreportlist.txt
  """
}
