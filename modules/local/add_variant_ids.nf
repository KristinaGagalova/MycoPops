/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ADD_VARIANT_IDS                                          (12.AddIDtoVCF.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Gives every SNP an ID of the form CHROM_POS_REF_ALT, e.g. ctg1_15000_C_A, so
    that downstream tools which track variants by ID (PLINK, LD pruning,
    ADMIXTURE) have a stable, unique handle on each site. PLINK in particular
    needs non-missing IDs.

        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT'

    The leading '+' matters: it fills in only the IDs that are currently missing
    and leaves any existing ID untouched (verified - an 'rs999' ID survives with
    '+' and is overwritten without it).

    This is the one step that needs bcftools rather than GATK, since GATK has no
    equivalent of --set-id.

    Output:
      <prefix>.snps_biallelic_withIDs.vcf.gz(.tbi)   FINAL variant set
      <prefix>.counts.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ADD_VARIANT_IDS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.19"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_1' :
        'quay.io/biocontainers/bcftools:1.19--h8b25389_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.withIDs.vcf.gz"), path("*.withIDs.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.counts.tsv")                                  , emit: counts
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: "--set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools annotate \\
        ${args} \\
        --threads ${task.cpus} \\
        -Oz \\
        -o ${prefix}.withIDs.vcf.gz \\
        ${vcf}

    tabix -p vcf ${prefix}.withIDs.vcf.gz

    {
        printf "12_variants_final\\t%s\\n" "\$(bcftools view -H ${prefix}.withIDs.vcf.gz | wc -l)"
        printf "12_variants_missing_id\\t%s\\n" "\$(bcftools query -f '%ID\\n' ${prefix}.withIDs.vcf.gz | grep -c '^\\.\$' || true)"
    } > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.withIDs.vcf.gz
    touch ${prefix}.withIDs.vcf.gz.tbi
    printf "12_variants_final\\t0\\n" > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
