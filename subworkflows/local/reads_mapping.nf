include { BWAMEM2_INDEX          } from '../../modules/nf-core/bwamem2/index/main'
include { ALIGN_MARKDUP_COVERAGE } from './align_markdup_coverage'

workflow READS_MAPPING {
    take:
    ch_reads      // [ val(meta), [ fastq_1, fastq_2 ] ]
    ch_fasta      // [ val(meta_ref), path(fasta) ]
    genome_size   // val

    main:
    // Build bwa-mem2 index only
    BWAMEM2_INDEX(ch_fasta)
    // Convert reference channels to VALUE channels so they're reused per sample
    def ch_fasta_val = ch_fasta.first()
    def ch_index_val = BWAMEM2_INDEX.out.index.first()

    ALIGN_MARKDUP_COVERAGE(
        ch_reads,
        ch_fasta_val,
        ch_index_val,
        genome_size
    )

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
    // versions = ALIGN_MARKDUP_COVERAGE.out.versions.mix(BWAMEM2_INDEX.out.versions)
}
