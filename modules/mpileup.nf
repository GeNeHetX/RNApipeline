process bcftools_mpileup{
	publishDir "${params.outputdir}/Mpileup_output/mpileup_file", mode: 'copy', pattern: '*.mpileup.gz'
	publishDir "${params.outputdir}/Mpileup_output/bcftools_file", mode: 'copy',  pattern: '*.bcftools.vcf.gz'

	input:
	tuple val(sample), path(bam), path(bai)
	path ref
	path regions_bed

	output:
	path "*.mpileup.gz"
	path "*.bcftools.vcf.gz"
	tuple val(sample), path("*.bcftools.vcf.gz"), emit:vc_file

	script:
	"""
	bcftools mpileup -f ${ref}"/ref.fa" ${bam} -R ${regions_bed} -d 500 >  "${sample}.mpileup"
	bcftools call  "${sample}.mpileup" -m -v -P ${params.bcftools_prior} -Ov -o "${sample}.bcftools.vcf"

	gzip "${sample}.mpileup"
	gzip "${sample}.bcftools.vcf"
	"""
}
