process bcftools_mpileup {
    cpus 1
    memory { params.mpileup_memory }

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path regions_bed

    output:
    tuple val(sample), path("${sample}.bcftools.vcf.gz"), emit: vc_file
    tuple val(sample), path("${sample}.depth.tsv"), emit: depths

    script:
    """
    set -euo pipefail
    bcftools mpileup -f ${ref}/ref.fa ${bam} -R ${regions_bed} \\
        -d ${params.bcftools_max_depth} -q ${params.variant_min_mapping_quality} \\
        -Q ${params.variant_min_base_quality} -a FORMAT/AD,FORMAT/DP -Ob -o ${sample}.pileup.bcf
    bcftools query -f '%CHROM\\t%POS[\\t%DP]\\n' ${sample}.pileup.bcf > ${sample}.depth.tsv
    bcftools call ${sample}.pileup.bcf -m -v -P ${params.bcftools_prior} -Ou | \\
        bcftools norm -f ${ref}/ref.fa -m -any -Oz -o ${sample}.bcftools.vcf.gz
    """
}
