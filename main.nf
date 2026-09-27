#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { POP_ANALYSIS_FLOW     } from './workflows/population-analysis'
include { FILTER_ONLY_FLOW      } from './workflows/filter-only'

//
// Default: the full pipeline (reads -> filtered cohort VCF)
//
workflow {
    POP_ANALYSIS_FLOW()
}

//
// Variant filtering only, on a cohort VCF that already exists:
//   nextflow run main.nf -entry FILTER --vcf <cohort.vcf.gz> --fasta <ref>
// See workflows/filter-only.nf
//
workflow FILTER {
    FILTER_ONLY_FLOW()
}
