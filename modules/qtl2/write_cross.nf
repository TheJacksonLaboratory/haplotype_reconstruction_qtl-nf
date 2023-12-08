process WRITE_CROSS {

  cpus 8
  memory 100.GB
  time '02:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/geno_probs", pattern: "*.RData", mode:'copy'

  input:
  file(cross_elements)

  output:
  path("*.RData"), emit: cross

  script:
  log.info "----- Make Control File and Cross Object -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/writeControlFile.R \
	${projectDir}/${params.covar} \
	${params.CCDOdataDir} \
	${projectDir}/${params.sample_folder}/qtl2genos
	
  """
}
