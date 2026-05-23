/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ALIGN_MARKDUP_COVERAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs the core population genetics alignment pipeline:
      1. BWA-MEM2 alignment + samtools sort  (nf-core module)
      2. Picard MarkDuplicates               (nf-core module)
      3. Coverage statistics                  (local module)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM           } from '../../modules/nf-core/bwamem2/mem/main'
include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main'
include { COVERAGE_STATS        } from '../../modules/local/coverage_stats'

workflow ALIGN_MARKDUP_COVERAGE {

    take:
    ch_reads      // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]
    ch_fasta      // channel: [ val(meta_ref), path(fasta) ]
    ch_index      // channel: [ val(meta_ref), path(bwa-mem2 index directory) ]
    ch_fai        // channel: [ val(meta_ref), path(fai) ]
    genome_size   // val: integer genome size in bp

    main:
    def ch_versions = channel.empty()

    //
    // STEP 1: Align paired-end reads with bwa-mem2 and sort
    //
    def sort_bam = true
    BWAMEM2_MEM(
        ch_reads,
        ch_index,
        ch_fasta,
        sort_bam
    )

    //
    // STEP 2: Mark duplicates with Picard
    //   Picard expects: tuple val(meta), path(bam)
    //                   tuple val(meta2), path(fasta), path(fai)
    //
    def ch_fasta_fai = ch_fasta.join(ch_fai.map { meta, fai -> [meta, fai] })
        .map { meta, fasta, fai -> [meta, fasta, fai] }

    PICARD_MARKDUPLICATES(
        BWAMEM2_MEM.out.bam,
        ch_fasta_fai
    )

    //
    // STEP 3: Calculate coverage statistics
    //
    COVERAGE_STATS(
        PICARD_MARKDUPLICATES.out.bam,
        genome_size
    )
    ch_versions = ch_versions.mix(COVERAGE_STATS.out.versions.first())

    emit:
    bam          = PICARD_MARKDUPLICATES.out.bam       // channel: [ val(meta), path(bam) ]
    metrics      = PICARD_MARKDUPLICATES.out.metrics   // channel: [ val(meta), path(metrics) ]
    coverage     = COVERAGE_STATS.out.coverage         // channel: [ val(meta), path(tsv) ]
    versions     = ch_versions                         // channel: [ path(versions.yml) ]
}
