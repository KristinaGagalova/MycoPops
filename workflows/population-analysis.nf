nextflow.enable.dsl = 2

include { READS_MAPPING     } from '../subworkflows/local/reads_mapping'
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { VARIANT_CALLING   } from '../subworkflows/local/variant_calling'

workflow POP_ANALYSIS_FLOW {
    main:
    if (!params.input) { error "Missing required parameter: --input" }
    if (!params.fasta) { error "Missing required parameter: --fasta" }
    if (params.min_coverage == null || !params.min_coverage.toString().isNumber() || (params.min_coverage as double) < 0) {
        error "Invalid --min_coverage '${params.min_coverage}': must be a number >= 0"
    }
    if (!params.ploidy.toString().isInteger() || (params.ploidy as int) < 1) {
        error "Invalid --ploidy '${params.ploidy}': must be a whole number >= 1"
    }

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
        params.genome_size,
        params.min_coverage
    )

    //
    // Reference files for GATK (.fai, .dict; decompressed if .gz)
    //
    PREPARE_REFERENCE(ch_fasta)

    //
    // Per-sample GVCFs with HaplotypeCaller (samples passing --min_coverage only)
    //
    VARIANT_CALLING(
        READS_MAPPING.out.bam,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai,
        PREPARE_REFERENCE.out.dict
    )

    emit:
    bam        = READS_MAPPING.out.bam        // [ meta, bam, bai ] samples passing --min_coverage
    bam_all    = READS_MAPPING.out.bam_all
    excluded   = READS_MAPPING.out.excluded
    metrics    = READS_MAPPING.out.metrics
    coverage   = READS_MAPPING.out.coverage
    genome_cov = READS_MAPPING.out.genome_cov
    gvcf       = VARIANT_CALLING.out.gvcf
    gvcf_tbi   = VARIANT_CALLING.out.tbi
    vcf        = VARIANT_CALLING.out.vcf          // joint-genotyped cohort VCF
    vcf_tbi    = VARIANT_CALLING.out.vcf_tbi

}
