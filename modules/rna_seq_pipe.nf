
nextflow.enable.dsl=2


process doOnlySTARnCount {

	publishDir "${params.outputdir}/FeautureCounts_output", mode: 'copy'

	input :
	path index
	tuple val(sample), file(fqFile)


	output:
	// path "${sample}StarOutAligned.sortedByCoord.out.bam"
	path "*StarOutLog.final.out"
	path "*_fastqc.{zip,html}"
	path "*exonscount.txt"
	path "*genecount.txt"

	when:
	starcount == true
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
	--outFilterMultimapNmax $params.outFilterMultimapNmax --winAnchorMultimapNmax $params.winAnchorMultimapNmax \
	--alignSJoverhangMin $params.alignSJoverhangMin\
	--alignSJDBoverhangMin $params.alignSJDBoverhangMin \
	--alignIntronMin $params.alignIntronMin \
	--outFilterMatchNminOverLread $params.outFilterMatchNminOverLread \
	--outFilterScoreMinOverLread $params.outFilterScoreMinOverLread\
	--outFilterMismatchNmax $params.outFilterMismatchNmax \
	--outFilterMismatchNoverLmax $params.outFilterMismatchNoverLmax



		featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $params.featureCountP -s $params.strand -O -o $sample'exonscount.txt' -f -t 'exon' -g 'exon_id' $sample'StarOutAligned.sortedByCoord.out.bam'

		featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $params.featureCountP -s $params.strand -O -o $sample'genecount.txt' -t 'exon' -g 'gene_id'  $sample'StarOutAligned.sortedByCoord.out.bam'




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
	--outFilterMultimapNmax $params.winAnchorMultimapNmax --winAnchorMultimapNmax $params.winAnchorMultimapNmax \
	--alignSJoverhangMin $params.alignSJoverhangMin\
	--alignSJDBoverhangMin $params.alignSJDBoverhangMin \
	--alignIntronMin $params.alignIntronMin \
	--outFilterMatchNminOverLread $params.outFilterMatchNminOverLread \
	--outFilterScoreMinOverLread $params.outFilterScoreMinOverLread\
	--outFilterMismatchNmax $params.outFilterMismatchNmax \
	--outFilterMismatchNoverLmax $params.outFilterMismatchNoverLmax


	"""


}

process FCounts {
	publishDir "${params.outputdir}/FeautureCounts_output", mode: 'copy'

	input:
	path allbams
	path index2



	output:
	path "exonscount_QC.txt"
	path "genecount_QC.txt"
	path "exonscount.txt"
	path "genecount.txt"


	"""

	featureCounts -T $task.cpus -F GTF -a  $index2/ref.gtf  -s $params.strand -O -o exonscount -f -t 'exon' -g 'exon_id' $allbams
	featureCounts -T $task.cpus -F GTF -a  $index2/ref.gtf  -s $params.strand -O -o genecount -t 'exon' -g 'gene_id'  $allbams


	awk 'NR>1' genecount > genecount.tab
	(head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d ) > genecount1.tab
	sort -d $index2/geneInfo.tab >geneinfos.tab
	awk 'NR>1' geneinfos.tab > geneinfos1.tab
	awk 'BEGIN { print "GeneidReference\tGeneName\tGeneType" } { print }' geneinfos1.tab > geneinfos2.tab
	paste -d geneinfos2.tab genecount1.tab > genecount_matrix.tab

	mv exonscount.summary exonscount_QC1.txt
	mv genecount.summary genecount_QC1.txt
	mv exonscount exonscount1.txt
	mv genecount_matrix.tab genecount1.txt


	awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' exonscount_QC1.txt > exonscount_QC.txt
	awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' genecount_QC1.txt > genecount_QC.txt
	awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' genecount1.txt | awk '{gsub("Chr", "Chromosome");print}'|awk '{gsub("Geneid", "GeneId");print}' > genecount.txt
	awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' exonscount1.txt | awk '{gsub("Chr", "Chromosome");print}'|awk '{gsub("Geneid", "ExonId");print}' > exonscount.txt


	"""
}

process multiqc {
	publishDir "${params.outputdir}/MULTIqc_output", mode: 'copy'
    input:

    path files

    output:
    file "*.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}
