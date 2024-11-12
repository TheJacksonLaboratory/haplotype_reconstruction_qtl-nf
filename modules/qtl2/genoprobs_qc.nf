process GENOPROBS_QC {

  cpus 1
  memory {300.GB * task.attempt}
  time {12.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.sample_folder}/geno_probs", pattern: "*.RData", mode:'copy'

  input:
  tuple file(cross)

  output:
  tuple file(cross), file("pr_36state.RData"), file("pr_8state.RData"), file("maxmarg_numeric.RData"), file("nxos.RData"), file("pr_errorlod.RData"), emit: genoprob_qc

  script:
  log.info "----- Performing Genotyping Quality Control and Measuring Crossovers -----"

  """
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/calcGenoProbs.R ${cross}
  """
}
