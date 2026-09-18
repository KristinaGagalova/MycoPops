/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPARE_REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Builds the reference files GATK needs next to the FASTA:
      1. gunzip the FASTA if it ends in .gz   (nf-core gunzip)
         GATK cannot read plain-gzipped FASTA
      2. <fasta>.fai                          (nf-core samtools/faidx)
      3. <fasta_basename>.dict                (nf-core gatk4/createsequencedictionary)

    All outputs are value channels, so they can be reused for every sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP                         } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'

workflow PREPARE_REFERENCE {

    take:
    ch_fasta   // channel: [ val(meta), path(fasta or fasta.gz) ]

    main:
    //
    // STEP 1: Decompress if needed
    //
    def ch_input = ch_fasta.branch { _meta, fasta ->
        gz:    fasta.name.endsWith('.gz')
        plain: true
    }

    GUNZIP(ch_input.gz)

    def ch_fasta_plain = GUNZIP.out.gunzip
        .mix(ch_input.plain)
        .first()

    //
    // STEP 2: FASTA index (.fai)
    //
    def get_sizes = false
    SAMTOOLS_FAIDX(
        ch_fasta_plain.map { meta, fasta -> [ meta, fasta, [] ] },
        get_sizes
    )

    //
    // STEP 3: Sequence dictionary (.dict)
    //
    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_plain)

    emit:
    fasta = ch_fasta_plain                                  // value channel: [ val(meta), path(fasta) ]
    fai   = SAMTOOLS_FAIDX.out.fai                          // value channel: [ val(meta), path(fai) ]
    dict  = GATK4_CREATESEQUENCEDICTIONARY.out.dict         // value channel: [ val(meta), path(dict) ]
}
