/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: COVERAGE_FILTER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Drops samples whose genome-wide mean depth (from GENOME_COVERAGE) is below
    `min_coverage`, so that downstream analysis only receives passing BAMs.

    A sample is excluded when:
      - mean_depth < min_coverage, or
      - mean_depth is missing / NA (e.g. no reads aligned)

    Reports written to ${params.outdir}/coverage/:
      - coverage_filter_report.tsv : every sample with PASS / FAIL
      - excluded_samples.tsv       : excluded samples only, with the reason
                                     (header-only when nothing was excluded)
    Each excluded sample is also logged as a warning.

    Set min_coverage to 0 to keep every sample that has a coverage value.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COVERAGE_FILTER {

    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ]
    ch_coverage   // channel: [ val(meta), path(genome_coverage.tsv) ]
    min_coverage  // val: minimum mean depth (number)

    main:
    def min_cov = min_coverage as double
    def outdir  = "${params.outdir}/coverage"

    //
    // Read mean_depth from each sample's GENOME_COVERAGE table
    //
    def ch_depth = ch_coverage.map { meta, tsv ->
        def rows  = tsv.splitCsv(header: true, sep: '\t')
        def value = rows ? rows[0].mean_depth?.toString() : null
        def depth = (value != null && value.isNumber()) ? value.toDouble() : null
        [ meta, depth ]
    }

    //
    // Pair each BAM with its depth and split into pass / fail
    //
    def ch_checked = ch_bam
        .join(ch_depth, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, bam, bai, depth ->
            def pass   = depth != null && depth >= min_cov
            def reason = depth == null ? 'coverage not available' :
                         pass          ? ''                       :
                                         "mean depth ${fmtDepth(depth)}x below minimum ${min_coverage}x"
            [ meta, bam, bai, depth, pass, reason ]
        }

    def ch_split = ch_checked.branch { _meta, _bam, _bai, _depth, pass, _reason ->
        pass: pass
        fail: true
    }

    //
    // Reports
    //
    ch_checked
        .map { meta, _bam, _bai, depth, pass, _reason ->
            "${meta.id}\t${fmtDepth(depth)}\t${min_coverage}\t${pass ? 'PASS' : 'FAIL'}\n"
        }
        .collectFile(
            name:     'coverage_filter_report.tsv',
            seed:     "sample\tmean_depth\tmin_coverage\tstatus\n",
            storeDir: outdir,
            sort:     true
        )

    ch_split.fail
        .map { meta, _bam, _bai, depth, _pass, reason ->
            "${meta.id}\t${fmtDepth(depth)}\t${min_coverage}\t${reason}\n"
        }
        .ifEmpty('')   // still write a header-only file when nothing is excluded
        .collectFile(
            name:     'excluded_samples.tsv',
            seed:     "sample\tmean_depth\tmin_coverage\treason\n",
            storeDir: outdir,
            sort:     true
        )

    ch_split.fail.subscribe { meta, _bam, _bai, _depth, _pass, reason ->
        log.warn "[COVERAGE_FILTER] Excluding sample '${meta.id}': ${reason}"
    }

    emit:
    bam      = ch_split.pass.map { meta, bam, bai, _depth, _pass, _reason -> [ meta, bam, bai ] }  // channel: [ val(meta), path(bam), path(bai) ]
    excluded = ch_split.fail.map { meta, _bam, _bai, depth, _pass, reason -> [ meta, depth, reason ] }  // channel: [ val(meta), val(depth), val(reason) ]
}

def fmtDepth(depth) {
    depth == null ? 'NA' : String.format('%.3f', depth as double)
}
