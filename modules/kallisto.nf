nextflow.enable.dsl=2


process KallistoPE {
	cpus { params.max_threads }
	memory { params.kallisto_memory }

	publishDir "${params.outputdir}/Kallisto_output", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:


	file "${Sample}/abundance.tsv"
	file "${Sample}/run_info.json"
	tuple val(Sample), path("${Sample}/abundance.tsv"), path("${Sample}/run_info.json"), emit: quant


	when:
	params.kallisto.toString().toBoolean()
	script:
	def read1 = []
	def read2 = []
	fastqFile.eachWithIndex { read, i ->
		if (i % 2 == 0) read1 << read
		else read2 << read
	}
	if (!read1 || read1.size() != read2.size())
		throw new IllegalArgumentException("Paired-end sample ${Sample} must have matching R1/R2 files")
	def prepareReads = ''
	def read1Arg = read1[0].toString()
	def read2Arg = read2[0].toString()
	if (read1.size() > 1) {
		prepareReads = "cat > read1_files.txt <<'READ1_FILES'\n${read1.join('\n')}\nREAD1_FILES\n" +
			"cat > read2_files.txt <<'READ2_FILES'\n${read2.join('\n')}\nREAD2_FILES\n" +
			"while IFS= read -r fq; do cat -- \"\$fq\"; done < read1_files.txt > \"${Sample}.kallisto_R1.fastq.gz\"\n" +
			"while IFS= read -r fq; do cat -- \"\$fq\"; done < read2_files.txt > \"${Sample}.kallisto_R2.fastq.gz\"\n"
		read1Arg = "${Sample}.kallisto_R1.fastq.gz"
		read2Arg = "${Sample}.kallisto_R2.fastq.gz"
	}
	"""
	set -euo pipefail
	for fq in $fastqFile; do
	    if [ ! -r "\$fq" ] || [ ! -s "\$fq" ]; then
	        echo "FASTQ missing, empty or unreadable on \$(hostname): \$fq" >&2
	        exit 1
	    fi
	done

	${prepareReads}

	kallisto quant -i $idx/kallisto_index -t ${task.cpus} $params.kallistoStrand \
		-o "${Sample}" "${read1Arg}" "${read2Arg}"
	"""



}

process Kallisto_single_end {
	cpus { params.max_threads }
	memory { params.kallisto_memory }
	publishDir "${params.outputdir}/Kallisto_output", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:
	tuple val(Sample), path("${Sample}/abundance.tsv"), path("${Sample}/run_info.json"), emit: quant

	when:
	params.kallisto.toString().toBoolean()
	script:
	def prepareReads = ''
	def readArg = fastqFile[0].toString()
	if (fastqFile.size() > 1) {
		prepareReads = "cat > read_files.txt <<'READ_FILES'\n${fastqFile.join('\n')}\nREAD_FILES\n" +
			"while IFS= read -r fq; do cat -- \"\$fq\"; done < read_files.txt > \"${Sample}.kallisto_single.fastq.gz\"\n"
		readArg = "${Sample}.kallisto_single.fastq.gz"
	}
	"""
	set -euo pipefail
	for fq in $fastqFile; do
	    if [ ! -r "\$fq" ] || [ ! -s "\$fq" ]; then
	        echo "FASTQ missing, empty or unreadable on \$(hostname): \$fq" >&2
	        exit 1
	    fi
	done

	${prepareReads}

	kallisto quant --single -l ${params.read_len} -s ${params.read_sd} $params.kallistoStrand \
		-i $idx/kallisto_index -t ${task.cpus} -o ${Sample} "${readArg}"
	"""



}
