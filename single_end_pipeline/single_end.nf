Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix )] }
.set { samples_ch }


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
	--genomeDir StarIndexed --genomeFastaFiles $/tmp/fa  --sjdbOverhang 100 --sjdbGTFfile $/tmp/gtf2 --genomeSAindexNbases 11
	samtools faidx $fa
	java -jar $params.picard CreateSequenceDictionary R= $fa O= genome.dict
	
	"""
}



process doSTAR {
	
	input :
	path index from index_ch
	set val(sample), file(fqFile)from samples_ch

	
	output: 
	file "${sample}StarOutAligned.sortedByCoord.out.bam" into bam_chcd
	file "*StarOutLog.final.out" into log_ch
	
	"""
	gunzip -c $fqFile >/tmp/fq
	STAR --genomeDir $index \
	--readFilesIn $/tmp/fq\
	--outFileNamePrefix $sample'StarOut' \
	--runThreadN 3 \
	--sjdbGTFfile $params.refGTF \
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
	file"genecount_matrix.tab"
	
	
	""" 
	gunzip -c $gtf >/tmp/annot
	featureCounts -T 3  -a $/tmp/annot -o exonscount -p -s 2 -f -t exon -g exon_id  $allbams
	featureCounts -t 'exon' -g 'gene_id' -a $/tmp/annot -T 3 -o genecount $allbams
	awk 'NR>1' genecount >genecount.tab
	(head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d)>genecount1.tab
	sort -d $index2/geneInfo.tab > geneinfos.tab
	awk 'NR>1' geneinfos.tab>geneinfos1.tab
	awk 'BEGIN { print "Geneid\tGenename\tGenetype" } { print }' geneinfos1.tab > geneinfos2.tab
	paste geneinfos2.tab  genecount1.tab  > genecount_matrix.tab
	
	"""