nextflow.enable.dsl=2

process gatk_vc {
	publishDir "${params.outputdir}/GATK4_output", mode: 'copy'
	input :
	file bamfile 
	path ref_data

	output: 
	file "*.vcf" 
	

	"""
	##marking duplicates 
	java -jar $params.picard  MarkDuplicates \
		I= $bamfile\
      		O=${bamfile.baseName}marked_duplicates.bam \
      		M=${bamfile.baseName}marked_dup_metrics.txt && \
      	java -jar $params.picard SortSam \
      		I=${bamfile.baseName}marked_duplicates.bam\
      		O=${bamfile.baseName}marked_duplicates_sorted.bam \
      		SORT_ORDER=coordinate
 
	##SplitNCigar 
	java -jar $params.gatk SplitNCigarReads \
      		-R $ref_data/ref.fa\
      		-I ${bamfile.baseName}marked_duplicates_sorted.bam\
      		-O ${bamfile.baseName}.split.bam 
	##adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		I= ${bamfile.baseName}.split.bam \
		O= ${bamfile.baseName}.split.RG.bam\
		RGLB=lib2 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=1 
	
	##BaseRecalibrator 
	java -jar $params.gatk BaseRecalibrator \
	  	-I ${bamfile.baseName}.split.RG.bam\
	  	-R $ref_data/ref.fa  \
	  	--known-sites $ref_data/knowns_variants.vcf \
	  	--use-jdk-inflater true \
	  	--use-jdk-deflater true \
	  	-O ${bamfile.baseName}.split.RG.recal.data.table
	##ApplyBQSR
	java -jar $params.gatk ApplyBQSR \
		-R $ref_data/ref.fa \
		-I ${bamfile.baseName}.split.RG.bam \
		--bqsr-recal-file ${bamfile.baseName}.split.RG.recal.data.table\
		-O ${bamfile.baseName}.abqsr.bam
	
	java -jar $params.gatk HaplotypeCaller  \
   	-R $ref_data/ref.fa \
   	-I ${bamfile.baseName}.abqsr.bam\
   	-O ${bamfile.SimpleName}.vcf \
   	-bamout ${bamfile.SimpleName}.bamout.bam

	"""
}



process Vep{

	publishDir "${params.outputdir}/VEP_output", mode: 'copy'
	
	input: 
	path vcf 
	
	output: 
	path "*_annot.tab"
	path "*_summary.txt"
	
	script: 
	"""
	vep -i $vcf -o ${vcf.baseName}_annot.tab --offline --cache --dir_cache /opt/vep/.vep --sift p,s,b --variant_class --nearest gene --gene_phenotype --hgvs --protein --symbol --canonical --mane  --var_synonyms --tab 


	"""
}


