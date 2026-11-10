nextflow.enable.dsl=2

process Vep {
	cpus 1
	memory { params.vep_memory }

	input:
	tuple val(sample), path(vcf)
	path ref_data
	val suffix

	output:
	tuple val(sample), path("*_annot.vcf.gz"), emit: annotations
	//path "*_annot.tab"

	script:
	def cache = "${ref_data}/VEP"
	"""
	set -euo pipefail
	vep --help | grep -E 'ensembl-vep[[:space:]]*:[[:space:]]*${params.ensembl_release}[.]'
	vep -i $vcf \
		--format vcf \
		-o ${sample}_${suffix}_annot.vcf \
		--vcf --no_stats --force_overwrite --transcript_version \
		--fasta $ref_data/ref.fa \
		--offline \
		--cache \
		--cache_version ${params.ensembl_release} \
		--dir_cache ${cache} \
		--sift b \
		--variant_class \
		--gene_phenotype \
		--hgvs \
		--protein \
		--symbol \
		--canonical \
		--mane \
		${params.vep_species}

	bgzip "${sample}_${suffix}_annot.vcf"
	"""
}
