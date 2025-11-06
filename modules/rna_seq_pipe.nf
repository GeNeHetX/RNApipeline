nextflow.enable.dsl=2

process doOnlySTARnCount {

	publishDir "${params.outputdir}/FeatureCounts_output", mode: 'copy'

	input :
	path index
	tuple val(sample), file(fqFile)


	output:
	// path "${sample}StarOutAligned.sortedByCoord.out.bam"
	path "*StarOutLog.final.out"
	path "*_fastqc.{zip,html}"
	path "*exonscount.txt"
	path "*genecount.txt"

	script:
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



	featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $params.featureCountP -s $params.FeatureCountStrand -O -o $sample'exonscount.txt' -f -t 'exon' -g 'exon_id' $sample'StarOutAligned.sortedByCoord.out.bam'

	featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $params.featureCountP -s $params.FeatureCountStrand -O -o $sample'genecount.txt' -t 'exon' -g 'gene_id'  $sample'StarOutAligned.sortedByCoord.out.bam'




	"""


}

process fastqc {
	publishDir "${params.outputdir}/Fastqc_output", mode: 'copy'

	input :
	tuple val(sample), file(fqFile)

	output:
	path '*_fastqc.{zip,html}'

	when:
	params.fastqc == true
	
	script:
	"""
	#echo $fqFile
	fastqc -t $task.cpus -q $fqFile
	"""
}

process samtools_index {
	
	input :
	tuple val(sample), path(bamFile)

	output:
	tuple val(sample), path(bamFile), path ("${sample}*.bai"), emit : align_files
	
	when:
	params.star == true
	
	script:
	"""
	samtools index -b ${bamFile} -o "${bamFile}.bai"
	"""
}

process doSTAR {
	publishDir "${params.outputdir}/UnmappedReads_output", mode: 'copy', pattern: "Unmapped_${sample}_R{1,2}.fastq.gz"
	publishDir "${params.outputdir}/StarLog_output", mode: 'copy', pattern: "*StarOutLog.final.out"

	input :
	path index
	tuple val(sample), file(fqFile)

	when:
	params.star == true
	
	output:
	path "${sample}StarOutAligned.sortedByCoord.out.bam"
	path "*StarOutLog.final.out"
	//path "*_fastqc.{zip,html}"
	path "Unmapped_${sample}_R{1,2}.fastq.gz", optional: true
	tuple val(sample), path("${sample}StarOutAligned.sortedByCoord.out.bam"), emit : bam4bai


	script:	
	"""

	# fastqc -t $task.cpus -q $fqFile

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
	--outFilterMismatchNoverLmax $params.outFilterMismatchNoverLmax \
	--outReadsUnmapped Fastx

	# Compress and rename Unmapped read fastqs
	if [ -s "${sample}StarOutUnmapped.out.mate1" ]; then
		gzip -c "${sample}StarOutUnmapped.out.mate1" > "Unmapped_${sample}_R1.fastq.gz"
	else
		echo "⚠️ No unmapped reads for mate1 (file empty or missing)"
	fi

	if [ -s "${sample}StarOutUnmapped.out.mate2" ]; then
		gzip -c "${sample}StarOutUnmapped.out.mate2" > "Unmapped_${sample}_R2.fastq.gz"
	else
		echo "ℹ️ Single-end mode or no unmapped mate2"
	fi

	"""
}

process FCounts {
	publishDir "${params.outputdir}/FeatureCounts_output", mode: 'copy'

	input:
	path bam
	path index
	tuple val(sample), file(fqFile)
	val featureCountP // -p if paired-end


	output:
	//path "*exonscount_QC.txt"
	//path "*genecount_QC.txt"
	path "*exonscount.txt"
	path "*genecount.txt"

	"""
	featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $featureCountP -s $params.FeatureCountStrand -O -o $sample'exonscount.txt' -f -t 'exon' -g 'exon_id' $bam

	featureCounts -T $task.cpus -F GTF -a  $index/ref.gtf $featureCountP -s $params.FeatureCountStrand -O -o $sample'genecount.txt' -t 'exon' -g 'gene_id' $bam
	
	"""

	// awk 'NR>1' genecount > genecount.tab
	// (head -n 1 genecount.tab && tail -n +2 genecount.tab | sort -d ) > genecount1.tab
	// sort -d $index/geneInfo.tab >geneinfos.tab
	// awk 'NR>1' geneinfos.tab > geneinfos1.tab
	// awk 'BEGIN { print "GeneidReference\tGeneName\tGeneType" } { print }' geneinfos1.tab > geneinfos2.tab
	// paste -d geneinfos2.tab genecount1.tab > genecount_matrix.tab

	// mv exonscount.summary exonscount_QC1.txt
	// mv genecount.summary genecount_QC1.txt
	// mv exonscount exonscount1.txt
	// mv genecount_matrix.tab genecount1.txt


	// awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' exonscount_QC1.txt > exonscount_QC.txt
	// awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' genecount_QC1.txt > genecount_QC.txt
	// awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' genecount1.txt | awk '{gsub("Chr", "Chromosome");print}'|awk '{gsub("Geneid", "GeneId");print}' > genecount.txt
	// awk '{gsub("StarOutAligned.sortedByCoord.out.bam", "");print}' exonscount1.txt | awk '{gsub("Chr", "Chromosome");print}'|awk '{gsub("Geneid", "ExonId");print}' > exonscount.txt
}

process multiqc {
	publishDir "${params.outputdir}/Multiqc_output", mode: 'copy'
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

process samtools_depth{

	publishDir "${params.outputdir}/SamtoolsDepth_output", mode: 'copy', pattern: "*_depth.txt"

	input:
	tuple val(sample), file(bamFile)
	path bedfile

	output:	
	path "*_depth.txt"

	when:
	params.samtools_depth == true

	script:
	if (params.no_multimapped == true)
		"""
		samtools depth -a -H -b $bedfile $bamFile -Q 255 -o ${sample}_q255_depth.txt
		"""
	else
		"""
		samtools depth -a -H -b $bedfile $bamFile -Q 0 -o ${sample}_q0_depth.txt
		"""
}