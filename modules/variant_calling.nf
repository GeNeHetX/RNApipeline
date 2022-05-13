nextflow.enable.dsl=2

process PicardSamsorted{
	

	input :
	file bamfile 
	output:
	file "*_sorted.bam" 
	
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
	
	
	  
	"""



}




process gatk_vc {
	publishDir "${params.outputdir}/VCF_output", mode: 'copy'
	input :
	path ref_data
	file sortedbam   
	
	
	output: 
	file "*.vcf" 
	
	
	"""
	
 
	##SplitNCigar 
	java -jar $params.gatk SplitNCigarReads \
      		-R $ref_data/ref.fa\
      		-I $sortedbam\
      		-O ${sortedbam.baseName}_split.bam 
	##adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		I= ${sortedbam.baseName}_split.bam \
		O= ${sortedbam.baseName}_split_RG.bam\
		RGLB=lib2 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=1 
	
	##BaseRecalibrator 
	java -jar $params.gatk BaseRecalibrator \
	  	-I ${sortedbam.baseName}_split_RG.bam\
	  	-R $ref_data/ref.fa  \
	  	--known-sites $ref_data/knowns_variants.vcf \
	  	--use-jdk-inflater true \
	  	--use-jdk-deflater true \
	  	-O ${sortedbam.baseName}_split_RG_recal_data.table
	##ApplyBQSR
	java -jar $params.gatk ApplyBQSR \
		-R $ref_data/ref.fa \
		-I ${sortedbam.baseName}_split_RG.bam \
		--bqsr-recal-file ${sortedbam.baseName}_split_RG_recal_data.table\
		-O ${sortedbam.baseName}_abqsr.bam
	##Running HaplotypeCaller 
	java -jar $params.gatk  HaplotypeCaller   -R $ref_data/ref.fa -I ${sortedbam.baseName}_abqsr.bam -O ${sortedbam.baseName}.vcf 
	
	rename sed -e 's/\StarOutAligned.sortedByCoord.outmarked_duplicates_sorted//' ${sortedbam.baseName}.vcf




	"""
}


process Vep{

	publishDir "${params.outputdir}/VEP_output", mode: 'copy'
	
	input: 
	path vcf 
	
	output: 
	path "*.vcf"
	
	script: 
	"""
	vep -i $vcf -o ${vcf.baseName}_annot.vcf --offline --cache --dir_cache /opt/vep/.vep --sift p,s,b --variant_class --nearest gene --gene_phenotype --hgvs --protein --symbol --canonical --mane  --var_synonyms --vcf 


	"""
}



