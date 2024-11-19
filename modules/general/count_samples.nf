process COUNT_SAMPLES {

  cpus 1
  memory 5.GB
  time '30minutes'
  
  input:
  tuple path(finalreport_files), val(project_id), path(covar_file), val(cross_type)

  output:
  tuple path(finalreport_files), val(project_id), path(covar_file), val(cross_type), file("line_count.txt"), emit: foo

  script:
  log.info "----- Counting Samples: Project ${project_id} -----"

  """
  wc -l ${covar_file} | awk '{print int((\$1+10)/20)*20}' > line_count.txt
  """
}