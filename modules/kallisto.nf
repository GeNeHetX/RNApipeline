nextflow.enable.dsl=2


process KallistoPE {

	publishDir "${params.outputdir}/kallistoOut", mode: 'copy'

	input:
	path idx
	tuple val(Sample), file(fastqFile)


	output:

	file "${Sample}_fr/abundance.h5"
	file "${Sample}_fr/abundance.tsv"
	file "${Sample}_fr/run_info.json"
	file "${Sample}_rf/abundance.h5"
	file "${Sample}_rf/abundance.tsv"
	file "${Sample}_rf/run_info.json"

	when:
	kallisto == true
	"""

	kallisto quant  -i $idx/kalliso_index -t ${task.cpus} --fr-stranded  -o "${Sample}_fr" ${fastqFile} && \
	kallisto quant  -i $idx/kalliso_index -t ${task.cpus} --rf-stranded  -o "${Sample}_rf" ${fastqFile}


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
