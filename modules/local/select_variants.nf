/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SELECT_VARIANTS                             (10.SelectVariantsSNP_INDELs.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Splits the FILTER-tagged VCF from VARIANT_FILTRATION into PASS SNPs and PASS
    INDELs, discarding everything the hard filters flagged:

        gatk SelectVariants --select-type-to-include SNP   --exclude-filtered true
        gatk SelectVariants --select-type-to-include INDEL --exclude-filtered true

    The original script then counted both with bcftools; the counts are written
    here instead, so no second container is needed.

    Only the SNPs continue down the pipeline (FILTER_SNPS -> ADD_VARIANT_IDS),
    which is what the original scripts do. The INDELs are published for
    inspection but not used further.

    Note: a site carrying both a SNP allele and an indel allele has GATK type
    MIXED and is selected by NEITHER call, so such sites are dropped here. That
    is the behaviour of the original script.

    Outputs:
      <prefix>.snps_pass.vcf.gz(.tbi)     PASS SNPs   -> continue
      <prefix>.indels_pass.vcf.gz(.tbi)   PASS INDELs -> side output
      <prefix>.counts.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SELECT_VARIANTS {
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
    tuple val(meta), path("*.snps_pass.vcf.gz")  , path("*.snps_pass.vcf.gz.tbi")  , emit: snps
    tuple val(meta), path("*.indels_pass.vcf.gz"), path("*.indels_pass.vcf.gz.tbi"), emit: indels
    tuple val(meta), path("*.counts.tsv")                                          , emit: counts
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[SELECT_VARIANTS] Available memory not known - defaulting to 3GB.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    #
    # PASS SNPs
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --select-type-to-include SNP \\
        --exclude-filtered true \\
        ${args} \\
        --tmp-dir . \\
        --output ${prefix}.snps_pass.vcf.gz

    #
    # PASS INDELs
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --select-type-to-include INDEL \\
        --exclude-filtered true \\
        ${args} \\
        --tmp-dir . \\
        --output ${prefix}.indels_pass.vcf.gz

    {
        printf "10_snps_pass\\t%s\\n"   "\$(gzip -cd ${prefix}.snps_pass.vcf.gz   | grep -vc '^#')"
        printf "10_indels_pass\\t%s\\n" "\$(gzip -cd ${prefix}.indels_pass.vcf.gz | grep -vc '^#')"
    } > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${prefix}.snps_pass ${prefix}.indels_pass; do
        echo "" | gzip > \$f.vcf.gz
        touch \$f.vcf.gz.tbi
    done
    printf "10_snps_pass\\t0\\n10_indels_pass\\t0\\n" > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """
}
