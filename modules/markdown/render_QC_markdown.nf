process QC_REPORT {

  cpus 1
  memory 50.GB
  time 2.hour

  container 'docker://sjwidmay/haplotype_reconstruction_qtl_nf:qc_markdown'

  publishDir "${params.sample_folder}/report", pattern: "*", mode:'copy'

  input:
  tuple file(cross), file(geno_probs), file(allele_probs), file(viterbi), file(crossovers), file(genotyping_erros), file(x_ints), file(y_ints), file(sample_qc_data)

  output:
  tuple file("QC_markdown.html"), file("QC_markdown.Rmd"), file("sample_QC_result.csv"), emit: qc_markdown

  script:
  log.info "----- Rendering Quality Control Report -----"

  """
  ls ${projectDir}/bin/scripts/markdown/QC_profile_template.Rmd
  cat ${projectDir}/bin/scripts/markdown/QC_profile_template.Rmd > QC_markdown_working.Rmd
  Rscript --vanilla ${projectDir}/bin/scripts/markdown/render_markdown.R QC_markdown_working.Rmd ${projectDir}/${params.sample_folder}
  mv QC_markdown_working.html QC_markdown.html
  mv QC_markdown_working.Rmd QC_markdown.Rmd
  """
}
