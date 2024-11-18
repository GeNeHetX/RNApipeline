nextflow.enable.dsl=2

process gatk_vc {
	publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	input :
	file bamfile 
	path ref_data

	output: 
	file "*.vcf" 
	
	when:
	params.variant_calling == true

	"""
	## 4)Sorting bams
	java -jar $params.picard SortSam \
		I= $bamfile\
      	O=${bamfile.baseName}_sorted.bam \
      	SORT_ORDER=coordinate

	## 3') Alignement Summary step (CollectAlignmentSummaryMetrics et CollectInsertSizeMetrics)
	## (TODO)

	## MergeBamAlignment (voir si ça améliore la qualité ou pas)
	## 5) Marking duplicates 
	java -jar $params.picard  MarkDuplicates \
		I= ${bamfile.baseName}_sorted.bam \
      	O=${bamfile.baseName}marked_duplicates_sorted.bam \
      	M=${bamfile.baseName}marked_dup_metrics.txt && \
 
	## 6) SplitNCigar 
	java -jar $params.gatk SplitNCigarReads \
      		-R $ref_data/ref.fa\
      		-I ${bamfile.baseName}marked_duplicates_sorted.bam\
      		-O ${bamfile.baseName}.split.bam 
	
	## 5') Adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		I= ${bamfile.baseName}.split.bam \
		O= ${bamfile.baseName}.split.RG.bam\
		RGLB=lib2 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=1 
	
	## 7) BaseRecalibrator + ApplyBQSR
	java -jar $params.gatk BaseRecalibrator \
	  	-I ${bamfile.baseName}.split.RG.bam\
	  	-R $ref_data/ref.fa  \
	  	--known-sites $ref_data/knowns_variants.vcf \
	  	--use-jdk-inflater true \
	  	--use-jdk-deflater true \
	  	-O ${bamfile.baseName}.split.RG.recal.data.table
	java -jar $params.gatk ApplyBQSR \
		-R $ref_data/ref.fa \
		-I ${bamfile.baseName}.split.RG.bam \
		--bqsr-recal-file ${bamfile.baseName}.split.RG.recal.data.table\
		-O ${bamfile.baseName}.abqsr.bam
	## AnalyzeCovariates (Evaluate and compare base quality score recalibration (BQSR) tables)
	## (Optional)

	## 8) BedToIntervalList (Bed réduit sur exons)
	## 9) INtervalListTools

	## 10) HaplotypeCaller
	java -jar $params.gatk HaplotypeCaller  \
   	-R $ref_data/ref.fa \
   	-I ${bamfile.baseName}.abqsr.bam\
   	-O ${bamfile.SimpleName}.vcf \
   	-bamout ${bamfile.SimpleName}.bamout.bam

	## 11) MergeVCFs (si parallelisation utile)
	## 12) Tabix (voir si c'est utile)
	## 13) VariantFiltration (fonction conditionnel)
	"""
}




process Vep{

	publishDir "${params.outputdir}/VEP_output", mode: 'copy'
	
	input: 
	path vcf 
	
	output: 
	path "*_annot.vcf"
	//path "*_annot.tab"
	
	when:
	params.vep == true
	
	script: 
	"""
	vep -i $vcf -o ${vcf.baseName}_annot.vcf --vcf --offline --cache --dir_cache ${params.vep_cache} --sift p,s,b --variant_class --gene_phenotype --hgvs --protein --symbol --canonical --mane
	"""
}


