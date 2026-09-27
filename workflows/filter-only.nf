nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: FILTER_ONLY_FLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs ONLY the variant-filtering step (steps 8 + 9) on a cohort VCF that
    already exists, so filter thresholds can be re-tried without touching the
    alignment or genotyping steps:

        nextflow run main.nf -entry FILTER \
            --vcf results/genotyping/cohort.genotyped.vcf.gz \
            --fasta /path/to/reference.fna.gz \
            --min_genotype_rate 0.8 \
            --outdir results_gr0.8 \
            -profile pawsey_setonix,singularity

    Use this instead of a full -resume run when:
      - you want to compare several thresholds (write each to its own --outdir)
      - the Nextflow `work` directory is gone (e.g. purged from /scratch), so
        -resume can no longer reuse the genotyping results
      - you want to filter a VCF produced somewhere else entirely

    The VCF index (.tbi) is used if it sits next to the VCF, and built if not.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_REFERENCE      } from '../subworkflows/local/prepare_reference'
include { GATK4_INDEXFEATUREFILE } from '../modules/nf-core/gatk4/indexfeaturefile/main'
include { VARIANT_FILTERING      } from '../subworkflows/local/variant_filtering'

workflow FILTER_ONLY_FLOW {
    main:
    if (!params.vcf)   { error "Missing required parameter: --vcf (a cohort VCF to filter)" }
    if (!params.fasta) { error "Missing required parameter: --fasta" }
    if (params.min_genotype_rate == null || !params.min_genotype_rate.toString().isNumber() ||
        (params.min_genotype_rate as double) < 0 || (params.min_genotype_rate as double) > 1) {
        error "Invalid --min_genotype_rate '${params.min_genotype_rate}': must be a number between 0 and 1"
    }

    def vcf_file = file(params.vcf, checkIfExists: true)
    def meta_vcf = [ id: vcf_file.simpleName ]

    //
    // Reference files GATK needs (.fai, .dict; decompressed if .gz)
    //
    def meta_ref = [ id: file(params.fasta).baseName ]
    def ch_fasta = channel.of([ meta_ref, file(params.fasta, checkIfExists: true) ])

    PREPARE_REFERENCE(ch_fasta)

    //
    // Reuse the .tbi next to the VCF when there is one, otherwise build it
    //
    def tbi_file = file("${params.vcf}.tbi")
    def ch_vcf = null

    if (tbi_file.exists()) {
        log.info "[FILTER] Using existing VCF index: ${tbi_file}"
        ch_vcf = channel.of([ meta_vcf, vcf_file, tbi_file ])
    }
    else {
        log.info "[FILTER] No index found at ${tbi_file} - building one"
        GATK4_INDEXFEATUREFILE(channel.of([ meta_vcf, vcf_file ]))
        ch_vcf = GATK4_INDEXFEATUREFILE.out.index.map { meta, idx -> [ meta, vcf_file, idx ] }
    }

    //
    // Hard filters + genotype-rate filter (same subworkflow the full run uses)
    //
    VARIANT_FILTERING(
        ch_vcf,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai,
        PREPARE_REFERENCE.out.dict,
        params.min_genotype_rate
    )

    emit:
    filtered_vcf   = VARIANT_FILTERING.out.vcf
    flagged_vcf    = VARIANT_FILTERING.out.flagged
    pass_vcf       = VARIANT_FILTERING.out.pass
    filter_summary = VARIANT_FILTERING.out.summary
}
