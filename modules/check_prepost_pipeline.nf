nextflow.enable.dsl=2

process Create_md5{
    publishDir "${params.outputdir}/Md5_output", mode: 'copy'

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
    publishDir "${params.outputdir}/Md5_output", mode: 'copy'

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

process AggregateSTAR {
    publishDir "${params.outputdir}/Aggregated_output", mode: 'copy'

    input:
    path feature_counts
    path star_logs
    path aggregate_script

    output:
    path 'StarFCGeneCount_s*.tsv.gz', emit: gene_matrix
    path 'StarFCExonCount_s*.tsv.gz', emit: exon_matrix
    path 'StarMetrics_s*.tsv', emit: star_metrics

    script:
    """
    set -euo pipefail
    mkdir -p FeatureCounts_output StarLog_output
    for f in ${feature_counts}; do
     ln -s "\$(realpath -- "\$f")" "FeatureCounts_output/\$(basename -- "\$f")"
    done
    for f in ${star_logs}; do
     ln -s "\$(realpath -- "\$f")" "StarLog_output/\$(basename -- "\$f")"
    done

    Rscript --vanilla ${aggregate_script} "\$PWD" ${task.cpus}
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
    path 'status_process.txt'

    script:
    def steps = [star: true, fastqc: true,
        fcounts: true, kallisto: report.run_kallisto,
        samtools: report.run_samtools, mpileup: report.run_variant_calling,
        vep_mpileup: report.run_variant_calling]
    steps.findAll { step, enabled -> enabled }.each { step, enabled ->
        def actual = completed.findAll { it[0] == step }.collect { it[1] }.sort()
        if (actual != report.samples.sort(false))
            error "Incomplete or duplicate ${step} results: expected ${report.samples}, got ${actual}"
    }
    if (completed.count { it[0] == 'multiqc' } != 1)
        error 'Missing MultiQC result'
    if (report.run_aggregate && completed.count { it[0] == 'aggregate' } != 1)
        error 'Missing STAR/FeatureCounts aggregation result'
    def json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(report))
    """
    cat > ${params.runNumber}_pipeline_summary.json <<'REPORT'
${json}
REPORT
    echo OK > status_process.txt
    """
}
