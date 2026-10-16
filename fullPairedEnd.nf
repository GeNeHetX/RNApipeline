nextflow.enable.dsl=2


include {doSTAR; FCounts; samtools_depth; multiqc; doOnlySTARnCount; samtools_index; fastqc} from './modules/rna_seq_pipe.nf'
include {gatk_vc; Vep as Vep_gatk; Vep as Vep_mpileup} from './modules/variant_calling.nf'
include {bcftools_mpileup} from './modules/mpileup.nf'
include {KallistoPE; Kallisto_single_end} from './modules/kallisto.nf'
include {buildref} from './modules/index.nf'
include {Create_md5; Verify_md5; Check_samples; Write_summary} from './modules/check_prepost_pipeline.nf'


// Fonction utilitaire pour formater dynamiquement la version du pipeline
def makeVersionTag(pipeline_version, genome_path) {
    if (!genome_path) {
        return "RNAv${pipeline_version.replaceAll(/^V/, '').replaceAll(/\\.0$/, '')}_customRef"
    }

    // Nettoyage de la version du pipeline -> "RNAv1.7.0" devient "RNAv1.7"
    def versionTag = "RNAv${pipeline_version.replaceAll(/^V/, '').replaceAll(/\\.0$/, '')}"

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
    if (reference_path == 'no_ref') {
        buildref(params.fasta_ref, params.GTF, params.cdna, params.known_vcf)
        reference = buildref.out
    } else {
        reference = file(reference_path, checkIfExists: true)
    }

    completed = Channel.empty()
    qc_files = Channel.empty()
    if (flag(params.fastqc)) {
        fastqc(samples_ch)
        completed = completed.mix(fastqc.out.reports.map { sample, files -> ['fastqc', sample] })
        qc_files = qc_files.mix(fastqc.out.reports.map { sample, files -> files })
    }
    if (flag(params.star)) {
        doSTAR(reference, samples_ch)
        completed = completed.mix(doSTAR.out.bam4bai.map { sample, bam -> ['star', sample] })
        qc_files = qc_files.mix(doSTAR.out[1])
    }
    if (flag(params.fcounts)) {
        FCounts(doSTAR.out.bam4bai, reference, flag(params.single_end) ? '' : '-p')
        completed = completed.mix(FCounts.out.counts.map { sample, exons, genes -> ['fcounts', sample] })
        qc_files = qc_files.mix(FCounts.out.summaries)
    }
    if (flag(params.kallisto)) {
        if (flag(params.single_end)) {
            Kallisto_single_end(reference, samples_ch)
            quant = Kallisto_single_end.out.quant
        } else {
            KallistoPE(reference, samples_ch)
            quant = KallistoPE.out.quant
        }
        completed = completed.mix(quant.map { sample, abundance, info -> ['kallisto', sample] })
    }
    if (flag(params.samtools_depth)) {
        samtools_depth(doSTAR.out.bam4bai, file(params.bed, checkIfExists: true))
        completed = completed.mix(samtools_depth.out.depths.map { sample, depth -> ['samtools', sample] })
    }
    if (flag(params.mpileup)) {
        samtools_index(doSTAR.out.bam4bai)
    }
    if (flag(params.gatk4)) {
        gatk_vc(doSTAR.out.bam4bai, reference)
        completed = completed.mix(gatk_vc.out.vc_file.map { sample, vcf -> ['gatk', sample] })
        if (flag(params.vep)) {
            Vep_gatk(gatk_vc.out.vc_file, reference, 'gatk4')
            completed = completed.mix(Vep_gatk.out.annotations.map { sample, vcf -> ['vep_gatk', sample] })
        }
    }
    if (flag(params.mpileup)) {
        bcftools_mpileup(samtools_index.out.align_files, reference, file(params.bed, checkIfExists: true))
        completed = completed.mix(bcftools_mpileup.out.vc_file.map { sample, vcf -> ['mpileup', sample] })
        if (flag(params.vep)) {
            Vep_mpileup(bcftools_mpileup.out.vc_file, reference, 'mpileup')
            completed = completed.mix(Vep_mpileup.out.annotations.map { sample, vcf -> ['vep_mpileup', sample] })
        }
    }
    if (flag(params.multiqc)) {
        multiqc(qc_files.flatten().collect())
        completed = completed.mix(multiqc.out[0].map { report -> ['multiqc', 'ALL'] })
    }

    emit:
    done = completed
}

workflow Main {
    if (!(params.runNumber ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Invalid runNumber'
    reference_path = params.ref ?: (params.reference_root ? "${params.reference_root}/ensembl_v${params.ensembl_release}_GRCh38" : null)
    if (!reference_path) error 'Provide --ref, or use -c configExamples/IAC.config'
    if (!flag(params.star) && (flag(params.fcounts) || flag(params.gatk4) || flag(params.mpileup) || flag(params.samtools_depth)))
        error 'FeatureCounts and variant/depth workflows require --star true'
    if (flag(params.multiqc) && !flag(params.star) && !flag(params.fastqc)) error 'MultiQC needs STAR or FastQC'
    if ((flag(params.mpileup) || flag(params.samtools_depth)) && !params.bed) error 'Mpileup/depth requires --bed'
    if (flag(params.vep) && !flag(params.gatk4) && !flag(params.mpileup)) error 'VEP needs a variant caller'
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
        pipeline   : "RNApipeline NF - V1.8.0",
        date       : date.toString(),
        genome     : reference_path,
        ensembl_release: params.ensembl_release,
        version    : makeVersionTag("V1.8.0", reference_path),
        sequenceID : params.runNumber,
        samples    : records.collect { it[0] },
        input_files: records.collect { [sample: it[0], fastqs: it[1].collect { f -> f.toString() }] },
        fastq_path : params.sampleInputDir,
        outputdir  : params.outputdir,
        scriptDir  : params.scriptDir,
        routine    : flag(params.routine),
        single_end : flag(params.single_end),
        run_star   : flag(params.star),
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
        limitBAMsortRAM : params.limitBAMsortRAM,
        run_multiqc   : flag(params.multiqc),
        run_fastqc : flag(params.fastqc),
        run_fcounts: flag(params.fcounts),
        featureCountStrand  : params.FeatureCountStrand,
        run_kallisto: flag(params.kallisto),
        kallistoStrand  : params.kallistoStrand,
        read_len : params.read_len,
        read_sd : params.read_sd,
        run_gatk   : flag(params.gatk4),
        gatk_version  : params.gatk,
        run_samtools : flag(params.samtools_depth),
        run_mpileup: flag(params.mpileup),
        region_bed : params.bed,
        pval: params.pval,
        run_vep    : flag(params.vep),
        containers : params.containers
      ]
    }
    Write_summary(summary_ch, Analysis.out.done.toList())
}

workflow {
    Main()
}
