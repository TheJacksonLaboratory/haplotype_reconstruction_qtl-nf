#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import workflow
if (params.workflow == "sample_qc_haplotype_reconstructions"){
  include {HR_QC} from './workflows/sample_qc_haplotype_reconstructions'
}
if (params.workflow == "QTL_Mapping"){
  include {QTL} from './workflows/QTL_Mapping'
}

// Conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "sample_qc_haplotype_reconstructions"){
    HR_QC()
    }
  if (params.workflow == "QTL_Mapping"){
    QTL()
    }
}
