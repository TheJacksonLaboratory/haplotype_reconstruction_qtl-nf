process WRITE_CROSS {

  cpus 1
  memory 100.GB
  time '02:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/qtl2genos", pattern: "*.csv", mode:'copy'
  publishDir "${params.sample_folder}/qtl2genos", pattern: "*.fst", mode:'copy'

  input:
  tuple file(genos), file(foundergenos), file(gmaps), file(pmaps)

  output:
  path("*.RData"), emit: cross

  script:
  log.info "----- Make Control File and Cross Object -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/writeControlFile.R \
	${genos} 
    ${foundergenos} \
    ${gmaps} \
    ${pmaps} \
    ${projectDir}/${params.covar}
	
  """
}
