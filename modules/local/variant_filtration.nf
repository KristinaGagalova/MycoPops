/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT_FILTRATION                                  (8.VariantFiltration.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Hard-filters the joint-genotyped cohort VCF on QUAL, QD, MQ, FS and the three
    rank-sum tests. Variants are only FLAGGED in the FILTER column; nothing is
    removed here. SELECT_VARIANTS (script 10) does the removing.

    Filter expressions are set in conf/modules.config (ext.args).

    Note: a variant whose annotation is absent is NOT flagged by a filter that
    tests it. This is standard GATK behaviour and matches the original script.
    For haploid isolates the three rank-sum annotations are often missing, so
    check the flagged_* counts in the summary to see what each filter caught.

    Outputs:
      <prefix>.flagged.vcf.gz(.tbi)  every variant, FILTER tagged
      <prefix>.counts.tsv            variant count + how many each filter flagged
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VARIANT_FILTRATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.2.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.2.0--py310hdfd78af_1' :
        'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.flagged.vcf.gz"), path("*.flagged.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.counts.tsv")                                  , emit: counts
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[VARIANT_FILTRATION] Available memory not known - defaulting to 3GB.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantFiltration \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        ${args} \\
        --tmp-dir . \\
        --output ${prefix}.flagged.vcf.gz

    {
        printf "08_samples\\t%s\\n"        "\$(gzip -cd ${prefix}.flagged.vcf.gz | awk -F'\\t' '/^#CHROM/{print NF-9; exit}')"
        printf "08_variants_input\\t%s\\n" "\$(gzip -cd ${prefix}.flagged.vcf.gz | grep -vc '^#')"
        gzip -cd ${prefix}.flagged.vcf.gz \\
            | awk -F'\\t' '!/^#/ { n = split(\$7, f, ";"); for (i = 1; i <= n; i++) if (f[i] != "PASS" && f[i] != ".") c[f[i]]++ }
                           END { for (k in c) printf "08_flagged_%s\\t%s\\n", k, c[k] }' \\
            | sort
    } > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.flagged.vcf.gz
    touch ${prefix}.flagged.vcf.gz.tbi
    printf "08_variants_input\\t0\\n" > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """
}
