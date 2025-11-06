nextflow.enable.dsl=2

process Create_md5{
    publishDir "${params.outputdir}/Md5_output", mode: 'copy'

	input:
	path(fastqDir)
    path(sampleCsv)

	output:
    file("md5Tab.tsv")
    file("status_md5.txt")

    script:
	"""
	python3 ${projectDir}/PrePostScripts/createCheckSum.py -md5 md5Tab.tsv -p $fastqDir -c $sampleCsv -e fastq.gz fastq
	"""
}

process Verify_md5{
    publishDir "${params.outputdir}/Md5_output", mode: 'copy'

	input:
    path(fastqDir)
    path(md5tab)

    output:
    path(md5tab)
    file("status_md5.txt")

	when:
	params.routine == false

    script:
	"""
	python  ${projectDir}/PrePostScripts/verifyCheckSum.py -md5 $md5tab -dir $fastqDir 

	"""
}

process Check_samples{
    publishDir "${params.outputdir}", mode: 'copy'

	input:
	path(sampleCsv)
    path(fastqDir)
    path(check_sample_script) // Path to the check_sample.py script

    output:
    file("sampleChecked.csv")
    file("status_samples.txt")

    script:
    """
    python3 ${check_sample_script} -csv ${sampleCsv}  -p ${fastqDir} -se ${params.single_end}
    """
    
}

process Check_process{
    publishDir "${params.outputdir}", mode: 'copy'
    input:
	file csvtab
    path(multiqc_files) // Collects all MultiQC output files
    path(check_process_script)

    output:
    file("status_process.txt")

    script:
    """
    python3  ${check_process_script} -path_res ${params.outputdir} -path_csv ${csvtab}
    """
}
