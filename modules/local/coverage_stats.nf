process COVERAGE_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.23.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23.1--h4c18ab8_0' :
        'biocontainers/samtools:1.23.1--h4c18ab8_0' }"

    input:
    tuple val(meta), path(bam)
    val(genome_size)

    output:
    tuple val(meta), path("*.coverage.tsv"), emit: coverage
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mapped_and_len=\$(samtools stats ${bam} \\
        | awk '
            /^SN[[:space:]]+reads mapped:/   {m=\$4}
            /^SN[[:space:]]+average length:/ {l=\$4}
            END { if (m>0 && l>0) print m, l; else print 0, 0 }
        ')

    mapped_reads=\$(echo "\$mapped_and_len" | awk '{print \$1}')
    avg_len=\$(echo "\$mapped_and_len" | awk '{print \$2}')

    coverage=\$(awk -v n="\$mapped_reads" -v L="\$avg_len" -v G="${genome_size}" \\
        'BEGIN{ if (G>0) printf "%.3f", (n*L)/G; else print "NA" }')

    printf "%s\\t%s\\t%s\\t%s\\n" "${prefix}" "\$mapped_reads" "${genome_size}" "\$coverage" \\
        > ${prefix}.coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
