/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT_FILTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Filtering of a jointly genotyped cohort VCF, ported from:

      8.VariantFiltration.sh   hard filters on QUAL, QD, MQ, FS and the three
                               rank-sum tests; variants are FLAGGED, not removed
      (bridge)                 keep only the variants tagged PASS
      9.FilterGenotypeRate.sh  keep variants called in at least
                               `min_genotype_rate` of samples

    All three run inside the single VARIANT_FILTRATION module, which uses one
    GATK container and writes bgzipped, tabix-indexed VCFs directly.

    Takes an indexed cohort VCF, so it can be driven either by VARIANT_CALLING
    or by the '-entry FILTER' workflow on a VCF that already exists.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANT_FILTRATION } from '../../modules/local/variant_filtration'

workflow VARIANT_FILTERING {

    take:
    ch_vcf             // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta           // value channel: [ val(meta), path(fasta) ]
    ch_fai             // value channel: [ val(meta), path(fai) ]
    ch_dict            // value channel: [ val(meta), path(dict) ]
    min_genotype_rate  // val: minimum fraction of samples with a genotype call

    main:
    VARIANT_FILTRATION(
        ch_vcf,
        ch_fasta,
        ch_fai,
        ch_dict,
        min_genotype_rate
    )

    emit:
    vcf     = VARIANT_FILTRATION.out.vcf       // channel: [ val(meta), path(vcf), path(tbi) ] FINAL filtered set
    flagged = VARIANT_FILTRATION.out.flagged   // channel: [ val(meta), path(vcf), path(tbi) ] all variants, FILTER tagged
    pass    = VARIANT_FILTRATION.out.pass      // channel: [ val(meta), path(vcf), path(tbi) ] PASS only
    summary = VARIANT_FILTRATION.out.summary   // channel: [ val(meta), path(tsv) ]
}
