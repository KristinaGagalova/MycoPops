#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pop-analysis-flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Population genetics alignment pipeline:
      1. BWA-MEM2 paired-end alignment + samtools sort
      2. Picard MarkDuplicates
      3. Coverage statistics

    Input: A CSV/TSV manifest with columns: sample, fastq_1, fastq_2

    GitHub : https://github.com/KristinaGagalova/MycoPops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_MARKDUP_COVERAGE } from '../subworkflows/local/align_markdup_coverage'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POP_ANALYSIS_FLOW {

    main:
    // Parse the input manifest (CSV or TSV)
    def ch_reads = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: params.manifest_sep ?: ',')
        .map { row ->
            def meta = [id: row.sample]
            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = file(row.fastq_2, checkIfExists: true)
            [meta, [fastq_1, fastq_2]]
        }

    // Reference genome — nf-core modules expect [ [id: 'genome'], path ]
    def meta_ref = [id: file(params.fasta).baseName]
    def ch_fasta = channel.of([meta_ref, file(params.fasta, checkIfExists: true)])
    def ch_index = channel.of([meta_ref, file(params.bwamem2_index, checkIfExists: true)])
    def ch_fai   = channel.of([meta_ref, file("${params.fasta}.fai", checkIfExists: true)])

    // Run the core pipeline
    ALIGN_MARKDUP_COVERAGE(
        ch_reads,
        ch_fasta,
        ch_index,
        ch_fai,
        params.genome_size
    )

    // Collect all per-sample coverage into a single summary file
    ALIGN_MARKDUP_COVERAGE.out.coverage
        .map { _meta, tsv -> tsv }
        .collectFile(
            name: 'coverage_summary.tsv',
            storeDir: "${params.outdir}/coverage",
            keepHeader: false,
            sort: true
        )

    emit:
    bam      = ALIGN_MARKDUP_COVERAGE.out.bam
    metrics  = ALIGN_MARKDUP_COVERAGE.out.metrics
    coverage = ALIGN_MARKDUP_COVERAGE.out.coverage
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY POINT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    POP_ANALYSIS_FLOW()
}
