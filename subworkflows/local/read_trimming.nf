/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: READ_TRIMMING                                   (fastp_array.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Adapter and quality trimming of the raw paired-end reads with fastp, before
    anything is aligned.

    fastp writes R1 and R2 as two SEPARATE files (--out1 / --out2), exactly as
    the original script does with -o / -O, so the pair stays usable by
    bwa-mem2. Verified on real data: both files come out with the same number of
    reads, the same read names in the same order, and no orphans - when either
    mate of a pair fails a filter, fastp drops BOTH mates rather than leaving a
    single read behind.

    Reports (.html, .json, .log) are published to ${params.outdir}/fastp/.
    The trimmed FASTQs are not published by default, since they are large and
    are passed straight to the aligner; see conf/containers.config to change it.

    This step is off by default. Turn it on with --trim_reads, and give the
    pipeline RAW reads in the sample sheet. If your sample sheet already points
    at reads that were trimmed beforehand, leave it off so they are not trimmed
    twice.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP } from '../../modules/nf-core/fastp/main'

workflow READ_TRIMMING {

    take:
    ch_reads   // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]

    main:
    //
    // fastp needs meta.single_end to pick its paired-end branch. It is added
    // here and removed again below, so the meta that flows downstream is
    // exactly the one that came in.
    //
    def discard_trimmed_pass = false
    def save_trimmed_fail    = false
    def save_merged          = false

    FASTP(
        ch_reads.map { meta, reads -> [ meta + [ single_end: false ], reads, [] ] },
        discard_trimmed_pass,
        save_trimmed_fail,
        save_merged
    )

    def ch_trimmed = FASTP.out.reads
        .map { meta, reads -> [ meta.findAll { k, _v -> k != 'single_end' }, reads ] }

    emit:
    reads     = ch_trimmed         // channel: [ val(meta), [ path(R1), path(R2) ] ]
    json      = FASTP.out.json     // channel: [ val(meta), path(json) ]
    html      = FASTP.out.html     // channel: [ val(meta), path(html) ]
    trim_log  = FASTP.out.log      // channel: [ val(meta), path(log) ]  ('log' is reserved in Nextflow)
}
