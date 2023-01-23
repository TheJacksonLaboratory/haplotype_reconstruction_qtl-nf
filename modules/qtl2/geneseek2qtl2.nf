process GS_TO_QTL2 {

  //cpus 1
  //memory 15.GB
  //time '00:30:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/qtl2genos", pattern: "*.csv", mode:'copy'

  input:
  tuple path(alleleCodes), path(FinalReport)

  output:
  file('*_geno*.csv'), emit: qtl2genos
  file('*_int.csv'), emit: qtl2ints

  script:
  log.info "----- Convert FinalReport File to R/qtl2 Genotypes and Intensities -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/geneseek2qtl2.R ${alleleCodes} ${params.sample_folder}/${FinalReport}
  """
}
