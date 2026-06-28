#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { POP_ANALYSIS_FLOW } from './workflows/population-analysis'

workflow {
    POP_ANALYSIS_FLOW()
}
