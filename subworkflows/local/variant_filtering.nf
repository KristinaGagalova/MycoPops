/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT_FILTERING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Turns the raw joint-genotyped cohort VCF into the final analysis-ready SNP
    set, following the WPM scripts.

    NOTE ON ORDER: the scripts do NOT run in numeric order. Their input/output
    paths show the real chain is 8 -> 10 -> 9 -> 11 -> 12, because script 9 reads
    'snps_pass_variantfiltration.vcf.gz', which script 10 produces. So the
    genotype-rate filter is applied to PASS SNPs, not to all PASS variants.

      VARIANT_FILTRATION   (8)  flag variants failing the hard filters
      SELECT_VARIANTS      (10) keep PASS SNPs and PASS INDELs, separately
      FILTER_SNPS          (9)  SNPs called in >= min_genotype_rate of samples
                           (11) then biallelic SNPs only
      ADD_VARIANT_IDS      (12) set IDs to CHROM_POS_REF_ALT

    Only the SNPs continue past SELECT_VARIANTS, as in the original scripts; the
    PASS INDELs are published as a side output.

    All per-step counts are merged into
    ${params.outdir}/variant_filtering/filtering_summary.tsv
    (metrics are prefixed with the script number so they sort into order).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANT_FILTRATION } from '../../modules/local/variant_filtration'
include { SELECT_VARIANTS    } from '../../modules/local/select_variants'
include { FILTER_SNPS        } from '../../modules/local/filter_snps'
include { ADD_VARIANT_IDS    } from '../../modules/local/add_variant_ids'

workflow VARIANT_FILTERING {

    take:
    ch_vcf             // channel: [ val(meta), path(vcf), path(tbi) ] joint-genotyped cohort VCF
    ch_fasta           // value channel: [ val(meta), path(fasta) ]
    ch_fai             // value channel: [ val(meta), path(fai) ]
    ch_dict            // value channel: [ val(meta), path(dict) ]
    min_genotype_rate  // val: minimum fraction of samples with a genotype call

    main:
    //
    // STEP 8: flag variants failing the hard filters
    //
    VARIANT_FILTRATION(
        ch_vcf,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    //
    // STEP 10: PASS SNPs and PASS INDELs
    //
    SELECT_VARIANTS(
        VARIANT_FILTRATION.out.vcf,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    //
    // STEPS 9 + 11: genotype rate, then biallelic (SNPs only)
    //
    FILTER_SNPS(
        SELECT_VARIANTS.out.snps,
        ch_fasta,
        ch_fai,
        ch_dict,
        min_genotype_rate
    )

    //
    // STEP 12: CHROM_POS_REF_ALT variant IDs
    //
    ADD_VARIANT_IDS(
        FILTER_SNPS.out.vcf
    )

    //
    // One summary table for the whole chain
    //
    VARIANT_FILTRATION.out.counts
        .mix(SELECT_VARIANTS.out.counts, FILTER_SNPS.out.counts, ADD_VARIANT_IDS.out.counts)
        .flatMap { _meta, tsv -> tsv.readLines().collect { line -> line + '\n' } }
        .collectFile(
            name:     'filtering_summary.tsv',
            seed:     "metric\tvalue\n",
            storeDir: "${params.outdir}/variant_filtering",
            sort:     { line -> line }     // sort by metric name, so 08_ .. 12_ read in order
        )

    emit:
    vcf            = ADD_VARIANT_IDS.out.vcf        // [ meta, vcf, tbi ] FINAL SNP set, with IDs
    biallelic      = FILTER_SNPS.out.vcf            // [ meta, vcf, tbi ] before IDs were added
    genotype_rate  = FILTER_SNPS.out.genotype_rate  // [ meta, vcf, tbi ] after the genotype-rate filter
    snps_pass      = SELECT_VARIANTS.out.snps       // [ meta, vcf, tbi ] PASS SNPs
    indels_pass    = SELECT_VARIANTS.out.indels     // [ meta, vcf, tbi ] PASS INDELs (side output)
    flagged        = VARIANT_FILTRATION.out.vcf     // [ meta, vcf, tbi ] all variants, FILTER tagged
}
