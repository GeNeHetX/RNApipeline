nextflow.enable.dsl=2

process Create_md5{

	input:
	path(fastqDir)
    path(sampleCsv)
    path(checksum_script)

	output:
    file("md5Tab.tsv")
    file("status_md5.txt")

    script:
	"""
	python3 ${checksum_script} -md5 md5Tab.tsv -p $fastqDir -c $sampleCsv -e fastq.gz fastq
	"""
}

process Verify_md5{

	input:
    path(fastqDir)
    path(md5tab)
    path(verify_script)

    output:
    path(md5tab)
    file("status_md5.txt")

	when:
	!params.routine.toString().toBoolean()

    script:
	"""
	python3 ${verify_script} -md5 $md5tab -dir $fastqDir

	"""
}

process Check_samples{

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

// Wait for actual process outputs; publishDir copies can finish asynchronously.
process Write_summary {
    publishDir "${params.outputdir}", mode: 'copy'

    input:
    val report
    val completed

    output:
    path "${params.runNumber}_pipeline_summary.json"

    script:
    def steps = [star: true, fastqc: true,
        fcounts: true, kallisto: report.run_kallisto,
        mpileup: report.run_variant_calling,
        vep_mpileup: report.run_variant_calling, metrics: true,
        somalier: report.run_somalier, hotspots: report.run_variant_calling]
    steps.findAll { step, enabled -> enabled }.each { step, enabled ->
        def actual = completed.findAll { it[0] == step }.collect { it[1] }.sort()
        if (actual != report.samples.sort(false))
            error "Incomplete or duplicate ${step} results: expected ${report.samples}, got ${actual}"
    }
    def json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(report))
    """
    cat > ${params.runNumber}_pipeline_summary.json <<'REPORT'
${json}
REPORT
    """
}
