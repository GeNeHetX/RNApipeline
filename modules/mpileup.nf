process remove_multimapped_reads{
	
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input: 
	tuple val(sample), path(bam), path(bai)
	
	output: 
	tuple val(sample), path("${sample}_STAR_filtered.bam"), path("${sample}_STAR_filtered.bam.bai"), emit:filtered_bam
	
	when:
	params.mpileup == true
	params.no_multimapped == true

	script:
	"""
	samtools view -h -q 255 -b $bam > ${sample}_STAR_filtered.bam
    samtools index -b ${sample}_STAR_filtered.bam -o ${sample}_STAR_filtered.bam.bai
	"""
}

process bcftools_mpileup{

	publishDir "${params.outputdir}/mpileup_output/mpileup_file", mode: 'copy', pattern: '*.mpileup.gz'
	publishDir "${params.outputdir}/mpileup_output/bcftools_file", mode: 'copy',  pattern: '*.bcftools.vcf.gz'
	
	container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'
	
	input: 
	tuple val(sample), path(bam), path(bai)
	path ref
	path regions_bed
	
	output: 
	path "*.mpileup.gz"
	path "*.bcftools.vcf.gz"
	tuple val(sample), path("*.bcftools.vcf.gz"), emit:vc_file
	
	when:
	params.mpileup == true

	script: 
	"""
	bcftools mpileup -f ${ref}"/ref.fa" ${bam} -R ${regions_bed} -d 500 >  "${sample}.mpileup"
	bcftools call  "${sample}.mpileup" -cv -p ${params.pval} -Ov -o "${sample}.bcftools.vcf"
	
	gzip "${sample}.mpileup"
	gzip "${sample}.bcftools.vcf"
	"""
}