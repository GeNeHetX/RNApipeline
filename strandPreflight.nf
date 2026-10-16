nextflow.enable.dsl=2

params.sampleList = null
params.sampleInputDir = null
params.single_end = false
params.ensembl_release = 107
params.ref = null
params.reference_root = null
params.outputdir = System.getenv('NF_RUN_RESULTS_DIR') ?: 'results'
params.salmon_container = 'combinelab/salmon@sha256:accc5581f515e10834d9127e67c4803e0cc713f7812b9699d54ad99b66e48742'

process SalmonIndex {
    container params.salmon_container

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index -t $transcriptome -i salmon_index -p $task.cpus
    """
}

process SalmonQuant {
    container params.salmon_container

    publishDir "${params.outputdir}/Salmon_output", mode: 'copy'

    input:
    path index
    tuple val(sample), path(fastqFiles)

    output:
    tuple val(sample), path("${sample}"), emit: results

    script:
    def mates1 = []
    def mates2 = []
    fastqFiles.eachWithIndex { read, i ->
        if (params.single_end.toString().toBoolean() || i % 2 == 0) mates1 << read
        else mates2 << read
    }
    def reads = params.single_end.toString().toBoolean()
        ? "-r ${mates1.join(' ')}"
        : "-1 ${mates1.join(' ')} -2 ${mates2.join(' ')}"
    """
    salmon quant -i $index -l A $reads -p $task.cpus -o "$sample"
    """
}

process SalmonSummary {
    publishDir "${params.outputdir}", mode: 'copy'

    input:
    path sampleDirs
    path summaryScript

    output:
    path 'strandedness_summary.tsv'

    script:
    """
    python3 $summaryScript $sampleDirs > strandedness_summary.tsv
    """
}

workflow {
    if (!params.sampleList) error 'Provide --sampleList'
    transcriptome_path = params.transcriptome ?: (params.ref ? "${params.ref}/cdna.fa" :
        (params.reference_root ? "${params.reference_root}/ensembl_v${params.ensembl_release}_GRCh38/cdna.fa" : null))
    if (!transcriptome_path) error 'Provide --transcriptome or configure --ref/reference_root'

    laneRows = file(params.sampleList, checkIfExists: true).readLines()
        .findAll { it.trim() && !it.trim().startsWith('#') }
        .collect { line ->
            def fields = line.trim().split(/\s+/)
            if (fields.size() != 4) error "Expected: sample read1 read2 directory; got: ${line}"
            def (sample, read1, read2, directory) = fields
            def base = (params.sampleInputDir ?: directory).replaceAll('/+$', '')
            if (params.single_end.toString().toBoolean() ? read2 != '-' : read2 == '-')
                error "Invalid read2 for ${sample}"
            def reads = params.single_end.toString().toBoolean() ? [read1] : [read1, read2]
            tuple(sample, reads.collect { file("${base}/${it}", checkIfExists: base.startsWith('s3://')) })
        }
    samples = laneRows.groupBy { it[0] }.collect { sample, lanes ->
        tuple(sample, lanes.collectMany { it[1] })
    }
    if (!samples || samples.any { !(it[0] ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) })
        error 'The sample sheet is empty or contains an invalid sample ID'
    fastqs = samples.collectMany { it[1] }.collect { it.toString() }
    if (fastqs.toSet().size() != fastqs.size()) error 'A FASTQ is assigned more than once in the sample sheet'

    SalmonIndex(file(transcriptome_path, checkIfExists: true))
    index = SalmonIndex.out.first()
    SalmonQuant(index, Channel.fromList(samples))
    sampleDirs = SalmonQuant.out.results.map { sample, directory -> directory }.collect()
    SalmonSummary(sampleDirs, file("${projectDir}/PrePostScripts/salmon_summary.py", checkIfExists: true))
}
