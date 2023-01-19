#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import workflow
if (params.workflow == "SampleQC_Haplotype_Reconstruction"){
  include {QC_HAP} from './workflows/SampleQC_Haplotype_Reconstruction'
}
if (params.workflow == "QTL_Mapping"){
  include {QTL} from './workflows/QTL_Mapping'
}

// Conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "SampleQC_Haplotype_Reconstruction"){
    QC_HAP()
    }
  if (params.workflow == "QTL_Mapping"){
    QTL()
    }
}
