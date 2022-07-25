




process buildref{
	
	input: 
	val fasta 
	val gtf 
	val cdna 
	val vcf 
	
	
	output:
	path "refdata" 
	
	
	"""
	mkdir refdata
	chmod +rwx refdata
	
	##downoalding the file 
	wget -nv -O ref.fa.gz $fasta 
	wget -nv -O ref.gtf.gz $gtf  
	wget -nv -O cdna.fa.gz $cdna 
	wget -nv -O knowns_variants.vcf.gz $vcf 
	
	gunzip ref.fa.gz
	gunzip cdna.fa.gz
	gunzip ref.gtf.gz
	gunzip knowns_variants.vcf.gz
	##create an  samtools index (needed for GATK4)
	samtools faidx ref.fa
	##kallisto index 
	kallisto index -i kalliso_index cdna.fa
	##you should specify the path to your picard.jar file
	java -jar $params.picard CreateSequenceDictionary R= ref.fa O= ref.dict
	java -jar $params.gatk IndexFeatureFile -I knowns_variants.vcf
	mv ref.gtf refdata/ref.gtf
	mv ref.fa.fai refdata/ref.fa.fai
	mv ref.dict refdata/ref.dict
	mv knowns_variants.vcf  refdata/knowns_variants.vcf 
	mv knowns_variants.vcf.idx refdata/knowns_variants.vcf.idx
	mv ref.fa refdata/ref.fa
	mv kalliso_index refdata/kalliso_index
	mv cdna.fa refdata/cdna.fa
	
	STAR --runThreadN 16 --runMode genomeGenerate --genomeDir refdata --genomeFastaFiles refdata/ref.fa  --sjdbOverhang 100 --sjdbGTFfile refdata/ref.gtf  --genomeSAindexNbases 11 
	
	"""

}



