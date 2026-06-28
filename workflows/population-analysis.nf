nextflow.enable.dsl = 2

include { READS_MAPPING } from '../subworkflows/local/reads_mapping'

workflow POP_ANALYSIS_FLOW {
    main:
    if (!params.input) { error "Missing required parameter: --input" }
    if (!params.fasta) { error "Missing required parameter: --fasta" }

    def ch_reads = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: params.manifest_sep ?: ',')
        .map { row ->
            def meta = [id: row.sample]
            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = file(row.fastq_2, checkIfExists: true)
            [meta, [fastq_1, fastq_2]]
        }

    def meta_ref = [id: file(params.fasta).baseName]
    def ch_fasta = channel.of([meta_ref, file(params.fasta, checkIfExists: true)])

    READS_MAPPING(
        ch_reads,
        ch_fasta,
        params.genome_size
    )

    emit:
    bam      = READS_MAPPING.out.bam
    metrics  = READS_MAPPING.out.metrics
    coverage = READS_MAPPING.out.coverage
}
