nextflow.enable.dsl=2

process gatk_vc {
	publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile) 
	path ref_data

	output: 
	tuple val(sample), file("*.vcf"), emit: vc_file
	
	when:
	params.gatk4 == true

	"""
	## SET TMP WORKING DIRECTORY

	## 4)Sorting bams
	java -jar $params.picard SortSam \
		-I $bamfile\
      	-O ${sample}_sorted.bam \
      	-SORT_ORDER coordinate \
		--TMP_DIR $params.tmpdir

	## 3') Alignement Summary step (CollectAlignmentSummaryMetrics et CollectInsertSizeMetrics)
	## (TODO)

	## 5') Adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		-I ${sample}_sorted.bam \
		-O ${sample}_sorted.RG.bam\
		-RGLB lib2 \
		-RGPL illumina \
		-RGPU unit1 \
		-RGSM 1 

	## MergeBamAlignment (voir si ça améliore la qualité ou pas)
	## 5) Marking duplicates 
	java -jar $params.picard  MarkDuplicates \
		-I ${sample}_sorted.RG.bam\
      	-O ${sample}_marked_duplicates_sorted.RG.bam \
      	-M ${sample}_marked_dup_metrics.txt && \
 
	## 6) SplitNCigar 
	java -Djava.io.tmpdir=$params.tmpdir -jar $params.gatk SplitNCigarReads \
      		-R $ref_data/ref.fa\
      		-I ${sample}_marked_duplicates_sorted.RG.bam\
      		-O ${sample}_marked_duplicates_sorted.RG.split.bam
	
	## 7) BaseRecalibrator + ApplyBQSR
	java -jar $params.gatk BaseRecalibrator \
	  	-I ${sample}_marked_duplicates_sorted.RG.split.bam\
	  	-R $ref_data/ref.fa  \
	  	--known-sites $ref_data/knowns_variants.vcf \
	  	--use-jdk-inflater true \
	  	--use-jdk-deflater true \
	  	-O ${sample}.split.RG.recal.data.table
	java -jar $params.gatk ApplyBQSR \
		-R $ref_data/ref.fa \
		-I ${sample}_marked_duplicates_sorted.RG.split.bam \
		--bqsr-recal-file ${sample}.split.RG.recal.data.table\
		-O ${sample}.abqsr.bam
	## AnalyzeCovariates (Evaluate and compare base quality score recalibration (BQSR) tables)
	## (Optional)

	## 8) BedToIntervalList (Bed réduit sur exons)
	java -jar $params.gatk BedToIntervalList \
    -I $ref_data/ref.ExonsOnly.bed \
    -O exons.interval_list \
    -R $ref_data/ref.fa \
    --SEQUENCE_DICTIONARY $ref_data/ref.dict
	
	## 9) INtervalListTools

	## 10) HaplotypeCaller
	java -jar $params.gatk HaplotypeCaller  \
   	-R $ref_data/ref.fa \
   	-I ${sample}.abqsr.bam\
   	-O ${sample}.vcf \
	-L exons.interval_list \
   	-bamout ${sample}.bamout.bam

	## 11) MergeVCFs (si parallelisation utile)
	## 12) Tabix (voir si c'est utile)
	## 13) VariantFiltration (fonction conditionnel)

	rm -rf ${params.tmpdir}/*
	rm -rf ${sample}_sorted.bam ${sample}_sorted.RG.bam ${sample}_marked_duplicates_sorted.RG.bam ${sample}_marked_duplicates_sorted.RG.split.bam ${sample}.abqsr.bam
	"""
}

process SortingSam {
	// publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)
	path ref_data

	output: 
	tuple val(sample), path("${sample}_sorted.bam"), emit : sorted_bam
 
	when:
	params.gatk4 == true

	"""
	## SET TMP WORKING DIRECTORY

	## 4)Sorting bams
	java -jar $params.picard SortSam \
		-I $bamfile\
      	-O ${sample}_sorted.bam \
      	-SORT_ORDER coordinate \
		--TMP_DIR $params.tmpdir
    """
}

process AddOrReplaceReadGroups {
	//publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)

	output: 
	tuple val(sample), path("${sample}_sorted.RG.bam"), emit : RG_bam 
	
	when:
	params.gatk4 == true

	"""
	## 3') Alignement Summary step (CollectAlignmentSummaryMetrics et CollectInsertSizeMetrics)
	## (TODO)

	## 5') Adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		-I ${bamfile} \
		-O ${sample}_sorted.RG.bam\
		-RGLB lib2 \
		-RGPL illumina \
		-RGPU unit1 \
		-RGSM 1 
    """
}

process MarkDuplicates {
	//publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)

	output: 
	tuple val(sample), path("${sample}_sorted.RG.MD.bam"), emit : MD_bam 
	file "${sample}metrics.MD.txt"
	
	when:
	params.gatk4 == true

	"""
	## MergeBamAlignment (voir si ça améliore la qualité ou pas)
	## 5) Marking duplicates 
	java -jar $params.picard  MarkDuplicates \
		-I ${bamfile}\
      	-O ${sample}_sorted.RG.MD.bam \
      	-M ${sample}metrics.MD.txt
	"""
}

process SplitNCigar {
	//publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)
	path ref_data

	output: 
	tuple val(sample), path("${sample}_sorted.RG.MD.split.bam"), emit : split_bam 
	
	when:
	params.gatk4 == true

	"""
	## 6) SplitNCigar 
	java -Xmx12g -jar $params.gatk SplitNCigarReads \
      		-R $ref_data/ref.fa\
      		-I ${bamfile}\
      		-O ${sample}_sorted.RG.MD.split.bam
	"""
}

process BaseRecalibrator {
	publishDir "${params.outputdir}/GATK4_output/recal_data_table", mode: 'copy', pattern: "*.recal.data.table"
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)
	path ref_data

	output: 
	tuple val(sample), path(bamfile), path("${sample}.recal.data.table"), emit: recal_files
	
	when:
	params.gatk4 == true

	"""
	## 7) BaseRecalibrator + ApplyBQSR
	java -jar $params.gatk BaseRecalibrator \
	  	-I ${bamfile}\
	  	-R $ref_data/ref.fa  \
	  	--known-sites $ref_data/knowns_variants.vcf \
	  	--use-jdk-inflater true \
	  	--use-jdk-deflater true \
	  	-O ${sample}.recal.data.table
	"""
}

process ApplyBQSR {
	//publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile), path(data_table)
	path ref_data

	output: 
	tuple val(sample), path("${sample}.abqsr.bam"), emit : abqsr_bam 
	
	when:
	params.gatk4 == true

	"""
	java -jar $params.gatk ApplyBQSR \
		-R $ref_data/ref.fa \
		-I $bamfile\
		--bqsr-recal-file $data_table\
		-O ${sample}.abqsr.bam
	## AnalyzeCovariates (Evaluate and compare base quality score recalibration (BQSR) tables)
	## (Optional)

	"""
}

process BedToIntervalList {
	//publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	path ref_data

	output: 
	file "exons.interval_list" 
	
	when:
	params.gatk4 == true

	"""
	## 8) BedToIntervalList (Bed réduit sur exons)
	java -jar $params.gatk BedToIntervalList \
    -I $ref_data/ref.ExonsOnly.bed \
    -O exons.interval_list \
    -R $ref_data/ref.fa \
    --SEQUENCE_DICTIONARY $ref_data/ref.dict

	## 9) INtervalListTools
	"""
}
	
process HaplotypeCaller {
	publishDir "${params.outputdir}/GATK4_output/", mode: 'copy', pattern: "*.vcf"
	container 'genehetx/genehetx-rnaseq:v1.6.0'
	
	input :
	tuple val(sample), path(bamfile)
	path ref_data
	file bedfile

	output: 
	file "${sample}.vcf" 
	
	when:
	params.gatk4 == true

	"""
	## 10) HaplotypeCaller
	java -jar $params.gatk HaplotypeCaller  \
   	-R $ref_data/ref.fa \
   	-I $bamfile \
   	-O ${sample}.vcf \
   	-bamout ${sample}.bamout.bam

	## 11) MergeVCFs (si parallelisation utile)
	## 12) Tabix (voir si c'est utile)
	## 13) VariantFiltration (fonction conditionnel)

	rm -rf ${params.tmpdir}/*
	"""
}

process Vep {
	publishDir "${params.outputdir}/VEP_output", mode: 'copy', pattern: "*_annot.vcf"
	container 'quay.io/biocontainers/ensembl-vep:113.2--pl5321h2a3209d_0'
	
	input: 
	tuple val(sample), path(vcf) 
	path ref_data
	val suffix
	
	output: 
	file "*_annot.vcf"
	//path "*_annot.tab"
	
	when:
	params.vep == true
	
	script: 
	"""
	vep -i $vcf \
		--format vcf \
		-o ${sample}_${suffix}_annot.vcf \
		--vcf \
		--fasta $ref_data/ref.fa \
		--offline \
		--cache \
		--dir_cache ${params.vep_cache} \
		--sift p,s,b \
		--variant_class \
		--gene_phenotype \
		--hgvs \
		--protein \
		--symbol \
		--canonical \
		--mane

	#gzip "${sample}_${suffix}_annot.vcf"
	"""
}
