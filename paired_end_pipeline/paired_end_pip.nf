Channel
    .fromFilePairs("/mnt/d/fastq/*_{1,2}.fastq.gz").set {sample_ch} 
annotation=file(params.refGTF)
genome= file(params.refFasta)
variant_file=file(params.vcf)
fastaindex= file (params.indexfasta)

process buildIndex {
	
	input:
	file(fasta) from genome
	file(gtf) from annotation
	
	output: 
	path "StarIndexed" into (index_ch, index_ch2)
	file "*.fai" into fai_idx_ch 
	file "*.dict" into dict_fa_ch 
	
	
	"""
	mkdir StarIndexed
	chmod +rwx StarIndexed
	gzip -c $fasta > /tmp/fa
	gzip -c $gtf > /tmp/gtf2
	STAR --runThreadN 3 --runMode genomeGenerate\
	--genomeDir StarIndexed --genomeFastaFiles /tmp/fa  --sjdbOverhang 100 --sjdbGTFfile /tmp/gtf2 --genomeSAindexNbases 11
	samtools faidx /tmp/fa
	java -jar $params.picard CreateSequenceDictionary R= /tmp/fa O= genome.dict
	
	"""
}



process doSTAR {
	
	input :
	path index from index_ch
	set val(sample1), file(fqFile)from samples_ch
	file gtf from annotation 

	
	output: 
	file "${sample1}StarOutAligned.sortedByCoord.out.bam" into(bam_chcd, bam_channel )
	file "*StarOutLog.final.out" into log_ch
	
	"""
	gunzip -c $fqFile >/tmp/fq1
	gunzip $gtf >/tmp/GTF
	STAR --genomeDir $index \
	--readFilesIn /tmp/fq1 \
	--outFileNamePrefix $sample1'StarOut' \
	--runThreadN 3 \
	--sjdbGTFfile /tmp/GTF \
	--twopassMode None --outFilterType BySJout  --seedSearchStartLmax 12 \
	--alignEndsType Local --outSAMtype BAM SortedByCoordinate \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--limitOutSJcollapsed 1000000 \
	--limitSjdbInsertNsj 1000000 \
	--outFilterMultimapNmax 100 --winAnchorMultimapNmax 50 \
	--alignSJoverhangMin 15 \
	--alignSJDBoverhangMin 1 \
	--alignIntronMin 20 \
	--outFilterMatchNminOverLread 0 \
	--outFilterScoreMinOverLread 0.3 \
	--outFilterMismatchNmax 33 \
	--outFilterMismatchNoverLmax 0.33 \
	
	"""

	
}

process FCounts {
	publishDir "${params.outputdir}/FC_output", mode: 'copy'

	input:
	file allbams from bam_chcd.collect()
	path index2 from index_ch2
	file gtf from annotation
	
	output: 
	file "exonscount.summary" into genecount_ch
	file"genecount.summary" into exoncount_ch
	file "genecount"
	file"exonscount" 
	file"genecount_matrix.tab" into gene_matrix_ch
	
	
	""" 
	gunzip -c $gtf >/tmp/annot
	featureCounts -T 3  -a /tmp/annot -o exonscount -p -s 2 -f -t exon -g exon_id  $allbams
	featureCounts -t 'exon' -g 'gene_id' -a /tmp/annot -T 3 -o genecount $allbams
	awk 'NR>1' genecount >genecount.tab
	(head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d)>genecount1.tab
	sort -d $index2/geneInfo.tab > geneinfos.tab
	awk 'NR>1' geneinfos.tab>geneinfos1.tab
	awk 'BEGIN { print "Geneid\tGenename\tGenetype" } { print }' geneinfos1.tab > geneinfos2.tab
	paste geneinfos2.tab  genecount1.tab  > genecount_matrix.tab
	
	"""

}





process PicardSamsorted{
	

	input :
	file bamfile from bam_channel
	output:
	file "*_sorted.bam" into samsorted
	when : 
	!params.no_gatk
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
	file sortedbam from samsorted 
	file variants from variant_file
	file fa from genome 
	file index fai_idx_ch 
	file dict from dict_fa_ch 
	
	output: 
	file "*.vcf" into vcf_ch
	file "*_out.bam" into bam_vcf_ch
	when:
	
	!params.no_gatk
	"""
	gunzip -c $fa >/tmp/fasta
	gunzip -c $variants >/tmp/vcf
	##SplitNCigar 
	java -jar $params.gatk SplitNCigarReads \
      -R /tmp/fasta \
      -I $sortedbam\
      -O ${sortedbam.baseName}_split.bam 
	##adding groups to reads 
	java -jar $params.picard AddOrReplaceReadGroups \
		I= ${sortedbam.baseName}_split.bam \
		O= ${sortedbam.baseName}_split_RG.bam\
		RGID=1 \
		RGLB=lib2 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=3 
	##BaseRecalibrator 
	java -jar $params.gatk BaseRecalibrator \
	  -I ${sortedbam.baseName}_split_RG.bam\
	  -R /tmp/fasta \
	  --known-sites $/tmp/vcf \
	  --use-jdk-inflater true \
	  --use-jdk-deflater true \
	  -O ${sortedbam.baseName}_split_RG_recal_data.table
	##ApplyBQSR
	java -jar $params.gatk ApplyBQSR \
		-R /tmp/fasta\
		-I ${sortedbam.baseName}_split_RG.bam \
		--bqsr-recal-file ${sortedbam.baseName}_split_RG_recal_data.table\
		-O ${sortedbam.baseName}_abqsr.bam
	##Running HaplotypeCaller 
	java -jar $params.gatk HaplotypeCaller -R /tmp/fasta -I ${sortedbam.baseName}_abqsr.bam -O ${sortedbam.baseName}.vcf -bamout ${sortedbam.baseName}_out.bam
	
	"""
}




