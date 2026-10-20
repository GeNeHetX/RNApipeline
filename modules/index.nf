process DownloadReference {
	input:
	val fasta
	val gtf
	val cdna
	val vcf
	val vep_cache_url
	val vep_cache_version
	val download_vep_cache

	output:
	tuple path("ref.fa"), path("ref.gtf"), path("cdna.fa"), path("knowns_variants.vcf"), path("VEP"), emit: files

	script:
	"""
	set -euo pipefail
	mkdir -p VEP

	wget --tries=3 --timeout=60 -nv -O ref.fa.gz "$fasta"
	wget --tries=3 --timeout=60 -nv -O ref.gtf.gz "$gtf"
	wget --tries=3 --timeout=60 -nv -O cdna.fa.gz "$cdna"
	wget --tries=3 --timeout=60 -nv -O knowns_variants.vcf.gz "$vcf"

	for file in ref.fa.gz ref.gtf.gz cdna.fa.gz knowns_variants.vcf.gz; do
		gzip -t "\$file"
	done
	gunzip -c ref.fa.gz > ref.fa
	gunzip -c ref.gtf.gz > ref.gtf
	gunzip -c cdna.fa.gz > cdna.fa
	gunzip -c knowns_variants.vcf.gz > knowns_variants.vcf
	rm -f ref.fa.gz ref.gtf.gz cdna.fa.gz knowns_variants.vcf.gz

	if [ "$download_vep_cache" = "true" ]; then
		wget --tries=3 --timeout=60 -nv -O vep_cache.tar.gz "$vep_cache_url"
		tar -xzf vep_cache.tar.gz -C VEP
		cache_dir=\$(find VEP -mindepth 2 -maxdepth 2 -type d -name '${vep_cache_version}_*' -print -quit)
		test -n "\$cache_dir"
		rm -f vep_cache.tar.gz
	fi

	test -s ref.fa
	test -s ref.gtf
	test -s cdna.fa
	test -s knowns_variants.vcf
	"""
}

process buildref {
	publishDir path: params.reference_build_root, mode: 'copy', overwrite: false, saveAs: { name ->
		name == 'refdata' ? (params.reference_name ?: "ensembl_v${params.ensembl_release}_GRCh38") : null
	}
	cpus { params.reference_build_threads }
	memory { params.reference_build_memory }

	input:
	tuple path(ref_fa), path(ref_gtf), path(cdna_fa), path(known_vcf), path(vep_cache)

	output:
	path "refdata"

	script:
	"""
	set -euo pipefail
	mkdir -p refdata

	# Build the small metadata indexes and the STAR/Kallisto indexes from the
	# downloaded files. GATK is used here only to index the known-variants VCF.
	samtools faidx "$ref_fa"
	kallisto index -i kallisto_index "$cdna_fa"
	java -jar $params.picard CreateSequenceDictionary R="$ref_fa" O=ref.dict
	java -jar $params.gatk IndexFeatureFile -I "$known_vcf"

	mv "$ref_fa" refdata/ref.fa
	mv "$ref_gtf" refdata/ref.gtf
	mv ref.fa.fai refdata/ref.fa.fai
	mv ref.dict refdata/ref.dict
	mv "$known_vcf" refdata/knowns_variants.vcf
	mv knowns_variants.vcf.idx refdata/knowns_variants.vcf.idx
	mv kallisto_index refdata/kallisto_index
	mv "$cdna_fa" refdata/cdna.fa
	mv "$vep_cache" refdata/VEP
		awk 'BEGIN { OFS="\\t" } \$3 == "exon" { print \$1, \$4 - 1, \$5 }' refdata/ref.gtf > refdata/ref.ExonsOnly.bed.raw
		awk 'BEGIN { OFS="\\t" } NR == 1 { chrom=\$1; start=\$2; end=\$3; next } \$1 == chrom && \$2 <= end { if (\$3 > end) end=\$3; next } { print chrom, start, end; chrom=\$1; start=\$2; end=\$3 } END { if (NR) print chrom, start, end }' refdata/ref.ExonsOnly.bed.raw > refdata/ref.ExonsOnly.bed
		rm refdata/ref.ExonsOnly.bed.raw

	STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
		--genomeDir refdata --genomeFastaFiles refdata/ref.fa \
		--sjdbOverhang 100 --sjdbGTFfile refdata/ref.gtf \
		--genomeSAindexNbases 11
	"""
}
