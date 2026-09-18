/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT_CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-sample calling and joint genotyping, ported from the WPM scripts:

      1. HaplotypeCaller per sample      (6.GATK_HaplotypeCaller_array.sh)
           -ERC GVCF -ploidy <ploidy>    -> <sample>.g.vcf.gz
      2. CombineGVCFs over all samples   (CombineGVCFs.sh)
                                         -> cohort.combined.g.vcf.gz
      3. IndexFeatureFile                (7.GenotypeGVCF.sh, index step)
                                         -> cohort.combined.g.vcf.gz.tbi
      4. GenotypeGVCFs                   (7.GenotypeGVCF.sh)
           --max-alternate-alleles 4     -> cohort.genotyped.vcf.gz

    Tool arguments are set in conf/modules.config.
    Only samples that passed COVERAGE_FILTER reach this subworkflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_HAPLOTYPECALLER  } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS     } from '../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_INDEXFEATUREFILE } from '../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_GENOTYPEGVCFS    } from '../../modules/nf-core/gatk4/genotypegvcfs/main'

workflow VARIANT_CALLING {

    take:
    ch_bam     // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta   // value channel: [ val(meta), path(fasta) ]
    ch_fai     // value channel: [ val(meta), path(fai) ]
    ch_dict    // value channel: [ val(meta), path(dict) ]

    main:
    //
    // STEP 1: HaplotypeCaller in GVCF mode, whole genome (no intervals, no DRAGstr model, no dbSNP)
    //
    GATK4_HAPLOTYPECALLER(
        ch_bam.map { meta, bam, bai -> [ meta, bam, bai, [], [] ] },
        ch_fasta,
        ch_fai,
        ch_dict,
        [ [:], [] ],   // dbsnp
        [ [:], [] ]    // dbsnp_tbi
    )

    //
    // STEP 2: Combine all per-sample GVCFs into one cohort GVCF.
    //   Sorted by sample id so the input order (and the -resume cache) is stable.
    //   Skipped entirely if no sample passed the coverage filter.
    //
    def ch_cohort = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .filter { samples -> samples.size() > 0 }
        .map { samples ->
            [
                [ id: 'cohort' ],
                samples.collect { s -> s[1] },   // GVCFs
                samples.collect { s -> s[2] }    // their .tbi indexes
            ]
        }

    GATK4_COMBINEGVCFS(
        ch_cohort,
        ch_fasta.map { _meta, fasta -> fasta },
        ch_fai.map   { _meta, fai   -> fai   },
        ch_dict.map  { _meta, dict  -> dict  }
    )

    //
    // STEP 3: Index the combined GVCF (the combinegvcfs module does not emit the .tbi)
    //
    GATK4_INDEXFEATUREFILE(
        GATK4_COMBINEGVCFS.out.combined_gvcf
    )

    //
    // STEP 4: Joint genotyping
    //
    GATK4_GENOTYPEGVCFS(
        GATK4_COMBINEGVCFS.out.combined_gvcf
            .join(GATK4_INDEXFEATUREFILE.out.index, failOnMismatch: true)
            .map { meta, gvcf, tbi -> [ meta, gvcf, tbi, [], [] ] },   // no intervals
        ch_fasta,
        ch_fai,
        ch_dict,
        [ [:], [] ],   // dbsnp
        [ [:], [] ]    // dbsnp_tbi
    )

    emit:
    gvcf          = GATK4_HAPLOTYPECALLER.out.vcf        // channel: [ val(meta), path(<sample>.g.vcf.gz) ]
    tbi           = GATK4_HAPLOTYPECALLER.out.tbi        // channel: [ val(meta), path(<sample>.g.vcf.gz.tbi) ]
    combined_gvcf = GATK4_COMBINEGVCFS.out.combined_gvcf // channel: [ val(meta), path(cohort.combined.g.vcf.gz) ]
    vcf           = GATK4_GENOTYPEGVCFS.out.vcf          // channel: [ val(meta), path(cohort.genotyped.vcf.gz) ]
    vcf_tbi       = GATK4_GENOTYPEGVCFS.out.tbi          // channel: [ val(meta), path(cohort.genotyped.vcf.gz.tbi) ]
}
