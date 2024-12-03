nextflow.enable.dsl=2


process KallistoPE {

	publishDir "${params.outputdir}/kallistoOut", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:


	file "${Sample}/abundance.tsv"
	file "${Sample}/run_info.json"


	when:
	params.kallisto == true
	"""


	kallisto quant  -i $idx/kallisto_index -t ${task.cpus} $params.kallistoStrand -o "${Sample}" ${fastqFile}


	"""



}

process Kallisto_single_end {
publishDir "${params.outputdir}/kallisto_output", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:
	file "${Sample}"

	when:
	kallisto == true
	"""

	kallisto quant --single -l ${params.read_len} -s ${params.read_sd}  -i $idx/kalliso_index -t ${task.cpus} -o ${Sample} ${fastqFile}
	"""



}
