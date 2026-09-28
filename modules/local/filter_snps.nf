/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_SNPS            (9.FilterGenotypeRate.sh + 11.biallelic_fitering.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two site-level filters on the PASS SNPs, in the order the original scripts
    run them:

      1. genotype rate   (9.FilterGenotypeRate.sh)
           keep SNPs called in at least `min_genotype_rate` of samples.
           Done with `gatk SelectVariants --max-nocall-fraction`, verified to
           give exactly the same sites as `vcftools --max-missing <rate>` at
           every threshold tested (1.0, 0.9, 0.8, 0.7, 0.5, 0.0), including a
           SNP missing in exactly 10% of samples.

      2. biallelic only  (11.biallelic_fitering.sh)
           keep SNPs with one REF and one ALT allele.
           Done with `gatk SelectVariants --restrict-alleles-to BIALLELIC`,
           verified to give the same sites as `vcftools --max-alleles 2` on this
           input. (The two differ only on monomorphic sites, which SELECT_VARIANTS
           has already removed.)

    Using GATK for both keeps this in one container and writes bgzipped,
    tabix-indexed VCFs directly, so the bgzip + tabix calls and the extra
    vcftools/htslib containers of the original scripts are not needed.

    Outputs:
      <prefix>.snps_genotyperate.vcf.gz(.tbi)   after filter 1
      <prefix>.snps_biallelic.vcf.gz(.tbi)      after filter 2
      <prefix>.counts.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_SNPS {
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
    tuple val(meta), path("*.snps_biallelic.vcf.gz")   , path("*.snps_biallelic.vcf.gz.tbi")   , emit: vcf
    tuple val(meta), path("*.snps_genotyperate.vcf.gz"), path("*.snps_genotyperate.vcf.gz.tbi"), emit: genotype_rate
    tuple val(meta), path("*.counts.tsv")                                                      , emit: counts
    path "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // vcftools --max-missing <rate>  ==  GATK --max-nocall-fraction (1 - rate)
    def rate       = min_genotype_rate as double
    def max_nocall = String.format('%.4f', 1.0d - rate)

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[FILTER_SNPS] Available memory not known - defaulting to 3GB.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    #
    # 1) genotype rate: keep SNPs called in >= ${rate} of samples
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --max-nocall-fraction ${max_nocall} \\
        --tmp-dir . \\
        --output ${prefix}.snps_genotyperate.vcf.gz

    #
    # 2) biallelic only
    #
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --reference ${fasta} \\
        --variant ${prefix}.snps_genotyperate.vcf.gz \\
        --restrict-alleles-to BIALLELIC \\
        ${args} \\
        --tmp-dir . \\
        --output ${prefix}.snps_biallelic.vcf.gz

    {
        printf "09_min_genotype_rate\\t${rate}\\n"
        printf "09_max_nocall_fraction\\t${max_nocall}\\n"
        printf "09_snps_genotyperate\\t%s\\n" "\$(gzip -cd ${prefix}.snps_genotyperate.vcf.gz | grep -vc '^#')"
        printf "11_snps_biallelic\\t%s\\n"    "\$(gzip -cd ${prefix}.snps_biallelic.vcf.gz    | grep -vc '^#')"
    } > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${prefix}.snps_genotyperate ${prefix}.snps_biallelic; do
        echo "" | gzip > \$f.vcf.gz
        touch \$f.vcf.gz.tbi
    done
    printf "09_snps_genotyperate\\t0\\n11_snps_biallelic\\t0\\n" > ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """
}
