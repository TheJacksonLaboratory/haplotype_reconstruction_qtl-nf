process QC_REPORT {

  cpus 1
  memory 50.GB
  time 2.hour

  container 'docker://sjwidmay/haplotype_reconstruction_qtl_nf:qc_markdown'

  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*", mode:'copy'

  input:
  tuple val(project_id), file(excluded_samples), file(genoprobs), file(alleleprobs), file(cross), file(viterbi), file(genotyping_errors), file(x_intensities), file(y_intensities), file(all_marker_intensities)

  output:
  tuple file("sample_QC.csv"), file("bad_markers.rds"), file("QC_markdown.html"), file("QC_markdown.Rmd"), emit: qc_markdown

  script:

  """
  ls ${projectDir}/bin/scripts/markdown/QC_profile_template.Rmd
  cat ${projectDir}/bin/scripts/markdown/QC_profile_template.Rmd > QC_markdown_working.Rmd
  Rscript --vanilla ${projectDir}/bin/scripts/markdown/render_markdown.R QC_markdown_working.Rmd
  mv QC_markdown_working.html QC_markdown.html
  mv QC_markdown_working.Rmd QC_markdown.Rmd
  """
}
