

process Kallisto_paired_end {
publishDir "${params.outputdir}/kallisto_output", mode: 'copy'

	input:
	path idx 
	tuple val(Sample), file(fastqFile)
	
	
	output: 
	
	file "${Sample}"
	"""
	
	kallisto quant -b ${params.bootstrap} -i $idx/kalliso_index -t ${task.cpus} -o ${Sample} ${fastqFile} 
	"""



}

process Kallisto_single_end {
publishDir "${params.outputdir}/kallisto_output", mode: 'copy'

	input:
	path idx 
	tuple val(Sample), file(fastqFile)
	
	
	output: 
	file "${Sample}"
	
	"""
	
	kallisto quant --single -l ${params.read_len} -s ${params.read_sd} -b ${params.bootstrap} -i $idx/kalliso_index -t ${task.cpus} -o ${Sample} ${fastqFile}
	"""



}

