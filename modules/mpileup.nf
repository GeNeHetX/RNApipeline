process bcftools_mpileup{

	publishDir "${params.outputdir}/mpileup_output/mpileup_file", mode: 'copy', pattern: '*.mpileup.gz'
	publishDir "${params.outputdir}/mpileup_output/bcftools_file", mode: 'copy',  pattern: '*.bcftools.vcf.gz'
	
	container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'
	
	input: 
	tuple val(sample), path(bam),path (bai)
	path ref
	path regions_bed
	
	output: 
	path "*.mpileup.gz"
	path "*.bcftools.vcf.gz"
	
	when:
	params.mpileup == true
	
	script: 
	"""
	bcftools mpileup -f ${ref}"/ref.fa" ${bam} -R ${regions_bed} -d 500 >  "${bam.baseName}.mpileup"
	bcftools call  "${bam.baseName}.mpileup" -cv -p ${params.pval} -Ov -o "${bam.baseName}.bcftools.vcf"
	
	gzip "${bam.baseName}.mpileup"
	gzip "${bam.baseName}.bcftools.vcf"
	"""
}