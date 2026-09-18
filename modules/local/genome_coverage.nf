/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME_COVERAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome-wide mean depth computed from aligned bases

    Unlike COVERAGE_STATS (mapped_reads * avg_read_length / genome_size), this:
      - needs no --genome_size: the reference length is taken from the BAM header
      - counts actual aligned bases, so soft-clipping and indels are handled
      - skips duplicates, secondary, QC-fail and unmapped reads
        (samtools depth default filter), so Picard-marked duplicates are excluded

    `-aa` reports every position of every contig, including contigs with no
    reads at all, so zero-coverage regions count towards the mean.

    Output columns (tab-separated, with header):
      sample  genome_length  aligned_bases  mean_depth  breadth_1x_pct
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENOME_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.23.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23.1--h4c18ab8_0' :
        'biocontainers/samtools:1.23.1--h4c18ab8_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.genome_coverage.tsv"), emit: coverage
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools depth \\
        -aa \\
        ${args} \\
        ${bam} \\
        | awk -v s="${prefix}" '
            BEGIN { OFS = "\\t" }
            { sum += \$3; if (\$3 > 0) cov++ }
            END {
                print "sample", "genome_length", "aligned_bases", "mean_depth", "breadth_1x_pct"
                if (NR > 0)
                    printf "%s\\t%d\\t%.0f\\t%.3f\\t%.2f\\n", s, NR, sum, sum / NR, 100 * cov / NR
                else
                    print s, 0, 0, "NA", "NA"
            }' \\
        > ${prefix}.genome_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample\\tgenome_length\\taligned_bases\\tmean_depth\\tbreadth_1x_pct\\n" > ${prefix}.genome_coverage.tsv
        printf "sample\\tgenome_length\\taligned_bases\\tmean_depth\\tbreadth_1x_pct\\n${prefix}\\t0\\t0\\t0.000\\t0.00\\n" > ${prefix}.genome_coverage.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
