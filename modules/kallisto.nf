nextflow.enable.dsl=2


process KallistoPE {

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
	"""


	for fq in $fastqFile; do
	    if [ ! -r "\$fq" ] || [ ! -s "\$fq" ]; then
	        echo "FASTQ missing, empty or unreadable on \$(hostname): \$fq" >&2
	        exit 1
	    fi
	done

	kallisto quant  -i $idx/kallisto_index -t ${task.cpus} $params.kallistoStrand -o "${Sample}" ${fastqFile}


	"""



}

process Kallisto_single_end {
	publishDir "${params.outputdir}/Kallisto_output", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:
	tuple val(Sample), path("${Sample}/abundance.tsv"), path("${Sample}/run_info.json"), emit: quant

	when:
	params.kallisto.toString().toBoolean()
	script:
	"""

	for fq in $fastqFile; do
	    if [ ! -r "\$fq" ] || [ ! -s "\$fq" ]; then
	        echo "FASTQ missing, empty or unreadable on \$(hostname): \$fq" >&2
	        exit 1
	    fi
	done

	kallisto quant --single -l ${params.read_len} -s ${params.read_sd} $params.kallistoStrand -i $idx/kallisto_index -t ${task.cpus} -o ${Sample} ${fastqFile}
	"""



}
