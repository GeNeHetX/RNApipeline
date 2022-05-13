
nextflow.enable.dsl=2





process buildIndex {
	publishDir "${params.outputdir}/staridx_output", mode: 'copy'
	
	input:
	path fasta 
	path gtf 
	path vcf 
	
	output: 
	path "ref_data" 
	
	
	"""
	
	
	mkdir ref_data
	chmod +rwx ref_data
	 
	gunzip -c $fasta > /tmp/ref.fa
	
	
	
	gunzip -c $gtf >/tmp/ref.gtf
	
	
	##create an  samtools index
	samtools faidx /tmp/ref.fa
	
	
	##create a dict 
	java -jar $params.picard CreateSequenceDictionary R= /tmp/ref.fa O= ref.dict
	
	##unzip the vcf file 
	
	gunzip -c $vcf>/tmp/knowns_variants.vcf
	
	
	##indexinf the vcf file 
	java -jar $params.gatk IndexFeatureFile --native-pair-hmm-threads $task.cpus -I /tmp/knowns_variants.vcf 
	
	STAR --runThreadN $task.cpus --runMode genomeGenerate\
	--genomeDir ref_data --genomeFastaFiles /tmp/ref.fa  --sjdbOverhang 100 --sjdbGTFfile /tmp/ref.gtf  --genomeSAindexNbases 11
	
	mv /tmp/ref.gtf ref_data/ref.gtf
	mv /tmp/ref.fa.fai ref_data/ref.fa.fai
	mv ref.dict ref_data/ref.dict
	mv /tmp/knowns_variants.vcf  ref_data/knowns_variants.vcf 
	mv /tmp/knowns_variants.vcf.idx ref_data/knowns_variants.vcf.idx
	mv /tmp/ref.fa ref_data/ref.fa
	
	"""
}



process doSTAR {
	
	input :
	
	
	path index 
	tuple val(sample), file(fqFile)
	

	
	output: 
	path "${sample}StarOutAligned.sortedByCoord.out.bam" 
	path "*StarOutLog.final.out" 
	path "*_fastqc.{zip,html}" 
	
	"""

	fastqc -t 8 -q $fqFile
	gunzip -c $fqFile|\
	STAR --genomeDir $index \
	--readFilesIn $fqFile\
	--readFilesCommand gunzip -c \
	--outFileNamePrefix $sample'StarOut' \
	--runThreadN $task.cpus \
	--sjdbGTFfile $index/ref.gtf\
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
	--outFilterMismatchNoverLmax 0.33 
	
	
	"""

	
}

process FCounts {
	publishDir "${params.outputdir}/FC_output", mode: 'copy'

	input:
	file allbams 
	path index2 

	
	
	output: 
	path "exonscount_QC.txt" 
	path "genecount_QC.txt" 
	path "exonscount.txt" 
	path "genecount_matrix.tab"
	
	
	""" 
	
	featureCounts -T $task.cpus -a  $index2/ref.gtf -M -o exonscount -f -t 'exon' -g 'exon_id' $allbams
	featureCounts -T $task.cpus -a  $index2/ref.gtf -M -o genecount -t 'exon' -g 'gene_id'  $allbams
	awk 'NR>1' genecount > genecount.tab 
	(head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d ) > genecount1.tab 
	sort -d $index2/geneInfo.tab >geneinfos.tab 
	awk 'NR>1' geneinfos.tab > geneinfos1.tab 
	awk 'BEGIN { print "Geneid_reference\tGene_name\tGene_type" } { print }' geneinfos1.tab > geneinfos2.tab
	paste -d geneinfos2.tab genecount1.tab > genecount_matrix.tab 
	
	mv exonscount.summary exonscount_QC.txt
	mv genecount.summary genecount_QC.txt
	mv exonscount exonscount.txt
	


	
	
	
	"""
}

process multiqc {
	publishDir "${params.outputdir}/MULTIqc_output", mode: 'copy'
    input:
   
    path files
    
    output:
    file "*.html"
    path "multiqc_data"
    path "QC_stats.txt"
    script:
    """
    multiqc .
   
    paste -d "\t" multiqc_data/multiqc_fastqc.txt  multiqc_data/multiqc_star.txt > QC_stats.txt
        """
}
