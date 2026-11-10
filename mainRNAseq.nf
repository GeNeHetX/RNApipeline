nextflow.enable.dsl=2


include {doSTAR; FCounts; samtools_index; fastqc} from './modules/rna_seq_pipe.nf'
include {Vep as Vep_mpileup} from './modules/variant_calling.nf'
include {bcftools_mpileup} from './modules/mpileup.nf'
include {SampleMetrics; Somalier; HotspotReport} from './modules/sample_outputs.nf'
include {KallistoPE; Kallisto_single_end} from './modules/kallisto.nf'
include {Create_md5; Verify_md5; Check_samples; Write_summary} from './modules/check_prepost_pipeline.nf'


// Fonction utilitaire pour formater dynamiquement la version du pipeline
def makeVersionTag(pipeline_version, genome_path) {
    if (!genome_path) {
        return "RNAv${pipeline_version.replaceAll(/^V/, '')}_customRef"
    }

    // Keep the full semantic pipeline version in the run metadata.
    def versionTag = "RNAv${pipeline_version.replaceAll(/^V/, '')}"

    // Extraction du nom de base du dossier du génome
    def genomeBase = genome_path.tokenize('/')[-1] ?: genome_path
    def genomeLower = genomeBase.toLowerCase()

    // Détection de la version Ensembl
    def m = (genomeBase =~ /ensembl_v(\d+)/)
    def genomeVersion = m.find() ? "Ensemblv${m.group(1)}" : "Custom"

    // Détection du type d’organisme
    def species = ""
    if (genomeLower.contains("humanmouse")) {
        species = "_PDX"
    } else if (genomeLower.contains("mouse")) {
        species = "_Mouse"
    } else {
        species = ""
    }

    // Construction finale
    return "${versionTag}_${genomeVersion}${species}"
}

// Nextflow 26 legacy CLI parameters can be strings, including "false".
def flag(value) { value.toString().toBoolean() }

// Keep each sample ID with all of its lane files. Never pair independent queues.
def validateSamples(rows) {
    if (!rows) error 'The sample sheet is empty'
    def ids = rows.collect { it[0] }
    if (ids.toSet().size() != ids.size()) error 'Duplicate sample IDs in the sample sheet'
    if (ids.any { !(it ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) }) error 'Invalid sample ID'
    def paths = rows.collectMany { it[1] }.collect { it.toString() }
    if (paths.toSet().size() != paths.size()) error 'A FASTQ is assigned more than once in the sample sheet'
    rows
}

workflow Analysis {
    take:
    samples_ch
    reference_path

    main:
    reference = file(reference_path, checkIfExists: true)

    completed = Channel.empty()
    fastqc(samples_ch)
    completed = completed.mix(fastqc.out.reports.map { sample, files -> ['fastqc', sample] })

    doSTAR(reference, samples_ch)
    completed = completed.mix(doSTAR.out.bam4bai.map { sample, bam -> ['star', sample] })
    FCounts(doSTAR.out.bam4bai, reference, flag(params.single_end) ? '' : '-p --countReadPairs')
    completed = completed.mix(FCounts.out.counts.map { sample, exons, genes -> ['fcounts', sample] })

    if (flag(params.kallisto)) {
        if (flag(params.single_end)) {
            Kallisto_single_end(reference, samples_ch)
            quant = Kallisto_single_end.out.quant
        } else {
            KallistoPE(reference, samples_ch)
            quant = KallistoPE.out.quant
        }
        completed = completed.mix(quant.map { sample, abundance, info -> ['kallisto', sample] })
        quant_info = quant.map { sample, abundance, info -> tuple(sample, info) }
    } else {
        quant_info = samples_ch.map { sample, reads -> tuple(sample, []) }
    }
    // One keyed bundle per sample, including an explicit disabled Kallisto record.
    metrics_inputs = samples_ch.map { sample, reads -> tuple(sample, reads.collect { it.name }) }
        .join(fastqc.out.reports, failOnDuplicate: true, failOnMismatch: true)
        .join(doSTAR.out.logs, failOnDuplicate: true, failOnMismatch: true)
        .join(FCounts.out.summaries, failOnDuplicate: true, failOnMismatch: true)
        .join(FCounts.out.counts.map { sample, exons, genes -> tuple(sample, genes) }, failOnDuplicate: true, failOnMismatch: true)
        .join(quant_info, failOnDuplicate: true, failOnMismatch: true)
    SampleMetrics(metrics_inputs, file("${params.scriptDir}/collect_metrics.py", checkIfExists: true))
    completed = completed.mix(SampleMetrics.out.metrics.map { sample, metrics -> ['metrics', sample] })

    if (flag(params.somalier) || flag(params.variant_calling)) {
        samtools_index(doSTAR.out.bam4bai)
    }
    if (flag(params.somalier)) {
        file("${reference_path}/somalier/sites.vcf.gz", checkIfExists: true)
        Somalier(samtools_index.out.align_files, reference)
        completed = completed.mix(Somalier.out.fingerprints.map { sample, fp -> ['somalier', sample] })
    }
    if (flag(params.variant_calling)) {
        file("${reference_path}/VEP", checkIfExists: true)
        file("${reference_path}/hotspots/validation.json", checkIfExists: true)
        bcftools_mpileup(samtools_index.out.align_files, reference,
            file("${reference_path}/hotspots/targets.bed", checkIfExists: true))
        completed = completed.mix(bcftools_mpileup.out.vc_file.map { sample, vcf -> ['mpileup', sample] })
        Vep_mpileup(bcftools_mpileup.out.vc_file, reference, 'mpileup')
        completed = completed.mix(Vep_mpileup.out.annotations.map { sample, vcf -> ['vep_mpileup', sample] })
        hotspot_inputs = Vep_mpileup.out.annotations.join(bcftools_mpileup.out.depths,
            failOnDuplicate: true, failOnMismatch: true)
        HotspotReport(hotspot_inputs, reference, file("${params.scriptDir}/hotspot_report.py", checkIfExists: true))
        completed = completed.mix(HotspotReport.out.reports.map { sample, report -> ['hotspots', sample] })
    }

    emit:
    done = completed
}

workflow Main {
    if (!(params.runNumber ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Invalid runNumber'
    reference_path = params.ref ?: (params.reference_root ? "${params.reference_root}/ensembl_v${params.ensembl_release}_GRCh38" : null)
    if (params.ref?.toString() == 'no_ref') error '`--ref no_ref` is not supported; build a reference with buildReference.nf, then provide its directory'
    if (!reference_path) error 'A prebuilt reference is required: set --ref or configure params.reference_root'
    [star: params.star, fastqc: params.fastqc, fcounts: params.fcounts].each { step, value ->
        if (value != null && !flag(value)) error "${step} is always on in mainRNAseq.nf; remove params.${step} = false"
    }
    if (params.mpileup != null || params.vep != null)
        error 'Use params.variant_calling to enable the complete bcftools and VEP workflow'
    if (!(params.bcftools_max_depth.toString() ==~ /[1-9]\d*/)) error '--bcftools_max_depth must be a positive integer'
    if (params.vep_cache_version && params.vep_cache_version.toString() != params.ensembl_release.toString())
        error 'VEP cache must use --ensembl_release; remove the separate --vep_cache_version override'
    reference_manifest = file("${reference_path}/reference_manifest.json")
    reference_info = reference_manifest.exists() ? new groovy.json.JsonSlurper().parseText(reference_manifest.text) : null
    if (reference_info && reference_info.ensembl_release.toString() != params.ensembl_release.toString())
        error 'Reference manifest release differs from --ensembl_release'
    if ((flag(params.somalier) || flag(params.variant_calling)) && !reference_info)
        error 'Somalier and hotspot calling require a 2.0 reference from buildReference.nf (reference_manifest.json missing)'
    if (params.sampleList && params.csvSample) error 'Provide either --sampleList or --csvSample'
    if (!params.sampleList && !params.csvSample) error 'Provide --sampleList or --csvSample'

    if (params.sampleList) {
        if (flag(params.check_md5)) error '--check_md5 is supported with legacy CSV input only'
        lane_rows = file(params.sampleList, checkIfExists: true).readLines()
            .findAll { it.trim() && !it.trim().startsWith('#') }
            .collect { line ->
                def fields = line.trim().split(/\s+/)
                if (fields.size() != 4) error "Expected: sample read1 read2 directory; got: ${line}"
                def (sample, read1, read2, directory) = fields
                def base = (params.sampleInputDir ?: directory).replaceAll('/+$', '')
                if (flag(params.single_end) ? read2 != '-' : read2 == '-') error "Invalid read2 for ${sample}"
                def reads = flag(params.single_end) ? [read1] : [read1, read2]
                tuple(sample, reads.collect { file("${base}/${it}", checkIfExists: base.startsWith('s3://')) })
            }
        rows = lane_rows.groupBy { it[0] }.collect { sample, lanes ->
            tuple(sample, lanes.collectMany { it[1] })
        }
        samples_ch = Channel.fromList(validateSamples(rows))
    } else {
        if (!params.sampleInputDir) error 'CSV input requires --sampleInputDir'
        if (flag(params.check_samples) && params.sampleInputDir.startsWith('s3://'))
            error 'S3 CSV input needs suffix columns and --check_samples false'
        if (flag(params.check_samples)) {
            Check_samples(file(params.csvSample), file(params.sampleInputDir),
                file("${params.scriptDir}/check_nb_sample.py"))
            sample_checked_csv = Check_samples.out[0]
        } else {
            sample_checked_csv = Channel.value(file(params.csvSample, checkIfExists: true))
        }
        samples_ch = sample_checked_csv.splitCsv(header: true).map { row ->
            if (!row.ID_Sample || row.suffix1 == null || (!flag(params.single_end) && row.suffix2 == null))
                error 'CSV needs ID_Sample,suffix1,suffix2; use --check_samples true to infer suffixes for local input'
            def suffixes = flag(params.single_end) ? [row.suffix1] : [row.suffix1, row.suffix2]
            tuple(row.ID_Sample, suffixes.collect {
                file("${params.sampleInputDir}/${row.ID_Sample}${it}.fastq.gz", checkIfExists: params.sampleInputDir.startsWith('s3://'))
            })
        }.toList().map { validateSamples(it) }.flatMap { it }
        if (flag(params.check_md5)) {
            if (flag(params.routine)) {
                Create_md5(file(params.sampleInputDir), sample_checked_csv,
                    file("${params.scriptDir}/createCheckSum.py"))
                checksum_status = Create_md5.out[1]
            } else {
                Verify_md5(file(params.sampleInputDir), file(params.md5),
                    file("${params.scriptDir}/verifyCheckSum.py"))
                checksum_status = Verify_md5.out[1]
            }
            verified = checksum_status.map { status ->
                if (status.text.trim() != 'OK') error 'FASTQ checksum verification failed'
                true
            }.first()
            samples_ch = samples_ch.combine(verified).map { sample, reads, ok -> tuple(sample, reads) }
        }
    }

    Analysis(samples_ch, reference_path)
    def date = new java.util.Date()
    summary_ch = samples_ch.toList().map { records ->
      [
        pipeline   : "RNApipeline NF - V${params.pipeline_version}",
        date       : date.toString(),
        genome     : reference_path,
        ensembl_release: params.ensembl_release,
        version    : makeVersionTag(params.pipeline_version, reference_path),
        sequenceID : params.runNumber,
        samples    : records.collect { it[0] },
        input_files: records.collect { [sample: it[0], fastqs: it[1].collect { f -> f.toString() }] },
        fastq_path : params.sampleInputDir,
        outputdir  : params.outputdir,
        scriptDir  : params.scriptDir,
        routine    : flag(params.routine),
        single_end : flag(params.single_end),
        run_star   : true,
        alignIntronMax : params.alignIntronMax,
        alignMatesGapMax  : params.alignMatesGapMax,
        limitOutSJcollapsed : params.limitOutSJcollapsed,
        limitSjdbInsertNsj  : params.limitSjdbInsertNsj,
        outFilterMultimapNmax : params.outFilterMultimapNmax,
        winAnchorMultimapNmax : params.winAnchorMultimapNmax,
        alignSJoverhangMin  : params.alignSJoverhangMin,
        alignSJDBoverhangMin  : params.alignSJDBoverhangMin,
        alignIntronMin  : params.alignIntronMin,
        outFilterMatchNminOverLread : params.outFilterMatchNminOverLread,
        outFilterScoreMinOverLread  : params.outFilterScoreMinOverLread,
        outFilterMismatchNmax : params.outFilterMismatchNmax,
        outFilterMismatchNoverLmax  : params.outFilterMismatchNoverLmax,
        run_fastqc : true,
        run_fcounts: true,
        featureCountStrand  : params.FeatureCountStrand,
        featureCountOptions : flag(params.single_end) ? '-O' : '-p --countReadPairs -O',
        run_kallisto: flag(params.kallisto),
        kallistoStrand  : params.kallistoStrand,
        fragment_length : params.read_len ?: params.fragment_length,
        read_sd : params.read_sd,
        run_somalier : flag(params.somalier),
        reference_manifest : reference_info,
        run_variant_calling: flag(params.variant_calling),
        run_mpileup: flag(params.variant_calling),
        bcftools_caller: '-m -v',
        region_bed : flag(params.variant_calling) ? "${reference_path}/hotspots/targets.bed" : null,
        bcftools_max_depth: params.bcftools_max_depth,
        variant_min_base_quality: params.variant_min_base_quality,
        variant_min_mapping_quality: params.variant_min_mapping_quality,
        bcftools_prior: params.bcftools_prior,
        run_vep    : flag(params.variant_calling),
        vep_cache_version: params.ensembl_release,
        containers : params.containers
      ]
    }
    Write_summary(summary_ch, Analysis.out.done.toList())
}

workflow {
    Main()
}
