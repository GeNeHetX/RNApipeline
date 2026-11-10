process SampleMetrics {
    publishDir "${params.outputdir}/Metrics_output", mode: 'copy'

    input:
    tuple val(sample), val(reads), path(qc), path(star_log), path(summaries), path(genes), path(kallisto_info)
    path metrics_script

    output:
    tuple val(sample), path("${sample}.metrics.tsv"), emit: metrics

    script:
    def kallisto_arg = params.kallisto.toString().toBoolean() ? "--kallisto ${kallisto_info}" : ''
    """
    python3 ${metrics_script} --sample ${sample} --reads ${reads.join(' ')} \
        --fastqc ${qc} --star ${star_log} --summaries ${summaries} --genes ${genes} \
        ${kallisto_arg} ${params.single_end.toString().toBoolean() ? '--single-end' : ''} \
        --output ${sample}.metrics.tsv
    """
}

process Somalier {
    publishDir "${params.outputdir}/Somalier_output", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.somalier"), emit: fingerprints

    script:
    """
    somalier extract --sites ${reference}/somalier/sites.vcf.gz -f ${reference}/ref.fa -d fingerprints ${bam}
    mv fingerprints/${sample}.somalier ${sample}.somalier
    """
}

process HotspotReport {
    publishDir "${params.outputdir}/Hotspot_output", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(depths)
    path reference
    path report_script

    output:
    tuple val(sample), path("${sample}.hotspots.tsv"), emit: reports

    script:
    """
    python3 ${report_script} --sample ${sample} --targets ${reference}/hotspots/targets.json \
        --depths ${depths} --vcf ${vcf} --output ${sample}.hotspots.tsv
    """
}
