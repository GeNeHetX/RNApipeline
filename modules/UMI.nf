nextflow.enable.dsl=2

process UMItools_extract {
	publishDir "${params.outputdir}/Fastq_process", mode: 'copy', pattern: "${sample}_process.{12}.fastq.gz"

	input:
	tuple val(sample), file(fqFile)

	output:	
	path val(sample), file(${sample}_process.2.fastq.gz), file(${sample}_process.1.fastq.gz)

	when:
	params.UMI == true

	script:
    """
	module load  umi_tools/1.1.2
    umi_tools extract --read2-in ${fqFile[0]} -I ${fqFile[1]} \
                  --bc-pattern=NNNNNNNN -L extract.log --stdout=${sample}_process.2.fastq.gz --read2-out=${sample}_process.1.fastq.gz
    """
}

process UMItools_dedup {
    //publishDir "${params.outputdir}/UMItools_stats", mode: 'copy', pattern: "${sample}_deduplicated_stats*"

	input:
	tuple val(sample), file(bamFile)

	output:	
	path val(sample), file(${sample}_deduplicated.bam)

	when:
	params.UMI == true

	script:
    """
    umi_tools dedup -I $bamFile --paired --output-stats=${sample}_deduplicated -S ${sample}_deduplicated.bam
    """
}