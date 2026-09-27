/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Everything variant-related, as two independent subworkflows:

      VARIANT_CALLING    per-sample GVCFs -> one jointly genotyped cohort VCF
                         (6.GATK_HaplotypeCaller_array.sh, CombineGVCFs.sh,
                          7.GenotypeGVCF.sh)

      VARIANT_FILTERING  cohort VCF -> hard-filtered, genotype-rate-filtered VCF
                         (8.VariantFiltration.sh, 9.FilterGenotypeRate.sh)

    They are kept separate so that filtering can be re-run on its own, without
    repeating the genotyping. Two ways to do that:

      - change only a filter setting and resubmit with -resume: Nextflow reuses
        the cached genotyping and re-runs the filtering task only
      - run the filtering on an existing cohort VCF, with no reads and no
        previous work directory:
            nextflow run main.nf -entry FILTER --vcf <cohort.vcf.gz> --fasta <ref>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANT_CALLING   } from './variant_calling'
include { VARIANT_FILTERING } from './variant_filtering'

workflow VARIANT {

    take:
    ch_bam             // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta           // value channel: [ val(meta), path(fasta) ]
    ch_fai             // value channel: [ val(meta), path(fai) ]
    ch_dict            // value channel: [ val(meta), path(dict) ]
    min_genotype_rate  // val: minimum fraction of samples with a genotype call

    main:
    //
    // Calling: reads -> cohort VCF
    //
    VARIANT_CALLING(
        ch_bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    //
    // Filtering: cohort VCF -> final filtered VCF
    //
    VARIANT_FILTERING(
        VARIANT_CALLING.out.vcf.join(VARIANT_CALLING.out.vcf_tbi, failOnMismatch: true),
        ch_fasta,
        ch_fai,
        ch_dict,
        min_genotype_rate
    )

    emit:
    // from VARIANT_CALLING
    gvcf           = VARIANT_CALLING.out.gvcf            // [ val(meta), path(<sample>.g.vcf.gz) ]
    gvcf_tbi       = VARIANT_CALLING.out.tbi             // [ val(meta), path(tbi) ]
    combined_gvcf  = VARIANT_CALLING.out.combined_gvcf   // [ val(meta), path(cohort.combined.g.vcf.gz) ]
    vcf            = VARIANT_CALLING.out.vcf             // [ val(meta), path(cohort.genotyped.vcf.gz) ] raw joint calls
    vcf_tbi        = VARIANT_CALLING.out.vcf_tbi         // [ val(meta), path(tbi) ]

    // from VARIANT_FILTERING
    filtered_vcf   = VARIANT_FILTERING.out.vcf           // [ val(meta), path(vcf), path(tbi) ] FINAL set
    flagged_vcf    = VARIANT_FILTERING.out.flagged       // [ val(meta), path(vcf), path(tbi) ] all variants, FILTER tagged
    pass_vcf       = VARIANT_FILTERING.out.pass          // [ val(meta), path(vcf), path(tbi) ] PASS only
    filter_summary = VARIANT_FILTERING.out.summary       // [ val(meta), path(tsv) ]
}
