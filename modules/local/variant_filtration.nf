/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT_FILTRATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Hard-filtering of the joint-genotyped cohort VCF. Ports two WPM scripts into
    one module, since both are single cohort-level jobs using the same reference:

      8.VariantFiltration.sh
        gatk VariantFiltration with 7 filter expressions (QUAL, QD, MQ, FS and
        the three rank-sum tests). This only FLAGS variants in the FILTER column;
        nothing is removed.            -> <prefix>.flagged.vcf.gz

      (bridge) gatk SelectVariants --exclude-filtered
        Keeps only variants tagged PASS. 9.FilterGenotypeRate.sh expects an
        already-PASS-filtered file as its input, so this step makes that explicit.
                                       -> <prefix>.pass.vcf.gz

      9.FilterGenotypeRate.sh
        Keeps variants genotyped in at least `min_genotype_rate` of samples.
        Done with `gatk SelectVariants --max-nocall-fraction`, which was verified
        to give exactly the same sites as `vcftools --max-missing <rate>` at every
        threshold tested (1.0, 0.9, 0.8, 0.7, 0.5, 0.0), including the boundary
        where a variant is missing in exactly 10% of samples.
        Using GATK keeps the whole module in one container and writes a bgzipped,
        tabix-indexed VCF directly, so the separate bgzip + tabix calls of the
        original script are not needed.
                                       -> <prefix>.filtered.vcf.gz   (FINAL)

    Also writes <prefix>.filter_summary.tsv: variant counts at each stage and how
    many variants each individual filter flagged.

    Notes:
      - A variant whose annotation is absent (e.g. no MQRankSum on a site with no
        heterozygous reads) is NOT flagged by a filter that tests it. This is
        standard GATK behaviour and matches the original script.
      - Filter expressions are set in conf/modules.config (ext.args), so they can
        be changed without editing this file.
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
    val min_genotype_rate

    output:
    tuple val(meta), path("*.flagged.vcf.gz")  , path("*.flagged.vcf.gz.tbi") , emit: flagged
    tuple val(meta), path("*.pass.vcf.gz")     , path("*.pass.vcf.gz.tbi")    , emit: pass
    tuple val(meta), path("*.filtered.vcf.gz") , path("*.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.filter_summary.tsv")                             , emit: summary
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''   // VariantFiltration: --filter-name / --filter-expression pairs
    def args2  = task.ext.args2 ?: ''   // extra SelectVariants options (e.g. --select-type-to-include SNP)
    def prefix = task.ext.prefix ?: "${meta.id}"

    // vcftools --max-missing <rate>  ==  GATK --max-nocall-fraction (1 - rate)
    def rate       = min_genotype_rate as double
    def max_nocall = String.format('%.4f', 1.0d - rate)

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[VARIANT_FILTRATION] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    #
    # STEP 1 (8.VariantFiltration.sh): flag variants failing the hard filters
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantFiltration \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        ${args} \\
        --tmp-dir . \\
        --output ${prefix}.flagged.vcf.gz

    #
    # STEP 2: keep only PASS variants
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${prefix}.flagged.vcf.gz \\
        --exclude-filtered \\
        ${args2} \\
        --tmp-dir . \\
        --output ${prefix}.pass.vcf.gz

    #
    # STEP 3 (9.FilterGenotypeRate.sh): keep variants called in >= ${rate} of samples
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${prefix}.pass.vcf.gz \\
        --max-nocall-fraction ${max_nocall} \\
        --tmp-dir . \\
        --output ${prefix}.filtered.vcf.gz

    #
    # Summary: variant counts per stage, and how many variants each filter flagged
    #
    {
        printf "metric\\tvalue\\n"
        printf "samples\\t%s\\n"                   "\$(gzip -cd ${prefix}.flagged.vcf.gz | awk -F'\\t' '/^#CHROM/{print NF-9; exit}')"
        printf "variants_input\\t%s\\n"            "\$(gzip -cd ${prefix}.flagged.vcf.gz  | grep -vc '^#')"
        printf "variants_pass_filters\\t%s\\n"     "\$(gzip -cd ${prefix}.pass.vcf.gz     | grep -vc '^#')"
        printf "variants_final\\t%s\\n"            "\$(gzip -cd ${prefix}.filtered.vcf.gz | grep -vc '^#')"
        printf "min_genotype_rate\\t${rate}\\n"
        printf "max_nocall_fraction\\t${max_nocall}\\n"
        gzip -cd ${prefix}.flagged.vcf.gz \\
            | awk -F'\\t' '!/^#/ { n = split(\$7, f, ";"); for (i = 1; i <= n; i++) if (f[i] != "PASS" && f[i] != ".") c[f[i]]++ }
                           END { for (k in c) printf "flagged_%s\\t%s\\n", k, c[k] }' \\
            | sort
    } > ${prefix}.filter_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${prefix}.flagged ${prefix}.pass ${prefix}.filtered; do
        echo "" | gzip > \$f.vcf.gz
        touch \$f.vcf.gz.tbi
    done
    printf "metric\\tvalue\\nvariants_final\\t0\\n" > ${prefix}.filter_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """
}
