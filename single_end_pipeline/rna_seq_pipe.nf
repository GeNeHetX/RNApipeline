
nextflow.enable.dsl=2

process doSTAR {
	
	input :
	
	
	path index 
	tuple val(sample), file(fqFile)

	
	output: 
	path "${sample}StarOutAligned.sortedByCoord.out.bam" 
	path "*StarOutLog.final.out" 
	path "*_fastqc.{zip,html}" 
	
	"""

	fastqc -t $task.cpus -q $fqFile
	
	STAR --genomeDir $index \
	--readFilesIn $fqFile\
	--readFilesCommand gunzip -c \
	--outFileNamePrefix $sample'StarOut' \
	--runThreadN $task.cpus \
	--sjdbGTFfile $index/ref.gtf\
	--twopassMode None --outFilterType BySJout  --seedSearchStartLmax 12 \
	--alignEndsType Local --outSAMtype BAM SortedByCoordinate \
	--alignIntronMax $params.alignIntronMax \
	--alignMatesGapMax $params.alignMatesGapMax \
	--limitOutSJcollapsed $params.limitOutSJcollapsed \
	--limitSjdbInsertNsj $params.limitSjdbInsertNsj \
	--outFilterMultimapNmax $params.outFilterMultimapNmax  --winAnchorMultimapNmax $params.winAnchorMultimapNmax \
	--alignSJoverhangMin params.alignSJoverhangMin  \
	--alignSJDBoverhangMin $params.alignSJDBoverhangMin \
	--alignIntronMin $params.alignIntronMin \
	--outFilterMatchNminOverLread $params.outFilterMatchNminOverLread \
	--outFilterScoreMinOverLread $params.outFilterScoreMinOverLread\
	--outFilterMismatchNmax $params.outFilterMismatchNmax \
	--outFilterMismatchNoverLmax $params.outFilterMismatchNoverLmax  
	
	
	"""

	
}

process FCounts {
	publishDir "${params.outputdir}/FC_output", mode: 'copy'

	input:
	file allbams 
	path index2 

	
	
	output: 
	path "exonscount.summary" 
	path "genecount.summary" 
	path "exonscount" 
	path "genecount_matrix.tab"
	
	
	""" 
	
	featureCounts -T $task.cpus -a  $index2/ref.gtf -M -o exonscount -f -t 'exon' -g 'exon_id' $allbams
	featureCounts -T $task.cpus -a  $index2/ref.gtf -M -o genecount -t 'exon' -g 'gene_id'  $allbams
	awk 'NR>1' genecount > genecount.tab 
	(head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d ) > genecount1.tab 
	sort -d $index2/geneInfo.tab >geneinfos.tab 
	awk 'NR>1' geneinfos.tab > geneinfos1.tab 
	awk 'BEGIN { print "Geneid_reference\tGenename\tGenetype" } { print }' geneinfos1.tab > geneinfos2.tab
	paste geneinfos2.tab genecount1.tab > genecount_matrix.tab 
	
	
	
	
	
	"""
}

process multiqc {
	publishDir "${params.outputdir}/MULTIqc_output", mode: 'copy'
    input:
   
    path 'fastqc/*'
    path logstar
    output:
    file "*.html"
    path "multiqc_data"
    path "QC_stats.txt"
    script:
    """
    multiqc .
   
    paste multiqc_data/multiqc_fastqc.txt  multiqc_data/multiqc_star.txt > QC_stats.txt

    """
}
