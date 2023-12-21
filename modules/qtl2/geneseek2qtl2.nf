process GS_TO_QTL2 {

  cpus 1
  memory {100.GB * task.attempt}
  time {5.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 3
  
  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/qtl2genos", pattern: "*.csv", mode:'copy'
  publishDir "${params.sample_folder}/qtl2genos", pattern: "*.fst", mode:'copy'

  input:
  path(FinalReport)

  output:
  path("*geno*.csv"), emit: qtl2genos
  path("*int.csv"), emit: qtl2ints
  path("*.fst"), emit: qtl2intsfst

  script:
  log.info "----- Convert FinalReport File to R/qtl2 Genotypes and Intensities -----"

  """
  echo ${FinalReport} > finalreportlist.txt

  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/geneseek2qtl2.R \
	${params.CCDOalleleCodes} \
	finalreportlist.txt \
	${projectDir}/${params.covar}
  """
}
