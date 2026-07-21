nextflow.enable.dsl=2


include {doSTAR; FCounts; samtools_depth; multiqc; doOnlySTARnCount; samtools_index; fastqc} from './modules/rna_seq_pipe.nf'
include {gatk_vc; Vep as Vep_gatk; Vep as Vep_deepvariant; Vep as Vep_mpileup} from './modules/variant_calling.nf'
include {Mosdepth; Bedtools; Deepvariant} from './modules/deepvariant.nf'
include {bcftools_mpileup} from './modules/mpileup.nf'
include {KallistoPE} from './modules/kallisto.nf'
include {buildref} from './modules/index.nf'
include {Create_md5; Verify_md5; Check_samples; Check_process; Write_summary} from './modules/check_prepost_pipeline.nf'
include {STAGE_WORKFLOW_RESOURCES} from './modules/workflow_resources.nf'


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

def referenceIdFromPath(reference_path) {
    if (!reference_path || reference_path.toString() == "no_ref") {
        return null
    }
    def normalized = reference_path.toString().replaceAll('/+$', '')
    return normalized.tokenize('/')[-1] ?: normalized
}

def referenceManifestFromPath(reference_path) {
    if (!reference_path || reference_path.toString() == "no_ref") {
        return null
    }
    def normalized = reference_path.toString().replaceAll('/+$', '')
    return "${normalized}/reference_manifest.json"
}

// Workflow pour l’analyse RNA-seq en paired-end
workflow Analysis_PE{
    take:

      samples_ch
      sample_csv // pas utilisé
      fastqDir // pas utilisé
      md5_error_file // sert à attendre les md5
      reference_path


    main:
      def featureCountP = " -p "  // paired end
      if(reference_path=="no_ref") {
        buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
        fastqc(samples_ch)
        doSTAR(buildref.out, samples_ch)
        FCounts(doSTAR.out[0], buildref.out, samples_ch, featureCountP)
        KallistoPE(buildref.out, samples_ch)
        gatk_vc(doSTAR.out.bam4bai, buildref.out)
      }
      else {
        fastqc(samples_ch)
        doSTAR(reference_path, samples_ch)
        KallistoPE(reference_path, samples_ch)

        if (params.mpileup || params.deepvariant) {
          samtools_index(doSTAR.out.bam4bai)
        }


        if (params.fcounts == true){
          FCounts(doSTAR.out[0],reference_path, samples_ch, featureCountP)
        }

        if (params.samtools_depth == true){
          samtools_depth(doSTAR.out.bam4bai, params.bed)
        }

        // VC with gatk4 + vep
        if (params.gatk4 == true){
          gatk_vc(doSTAR.out.bam4bai, reference_path)
          Vep_gatk(gatk_vc.out.vc_file, reference_path,"gatk4")
        }

        // VC with mpileup + vep
        if (params.mpileup == true){
          bcftools_mpileup(samtools_index.out.align_files , reference_path, params.bed)
          Vep_mpileup(bcftools_mpileup.out.vc_file, reference_path,"mpileup")
        }

        // VC with deepvariant + vep
        if (params.deepvariant == true){
          Mosdepth(doSTAR.out[0],samtools_index.out[0],samples_ch)
          Bedtools(Mosdepth.out[0],samples_ch)
          Deepvariant(doSTAR.out[0], samtools_index.out[0], samples_ch, Bedtools.out[0], reference_path, params.modelckptdeepar)
          Vep_deepvariant(Deepvariant.out, reference_path,"deepvariant")
        }


      }
      //Agregate quality results
      multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())

    emit:
      samples_ch
      qc_out = multiqc.out[1]

}

// Workflow pour l’analyse RNA-seq en single-end
workflow Analysis_SE{
      take:

      samples_ch
      md5_error_file // sert à attendre les md5
      reference_path


    main:
      def featureCountP = " "  // single end

      if(reference_path=="no_ref") {
        buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
        fastqc(samples_ch)
        doSTAR(buildref.out, samples_ch)
        FCounts(doSTAR.out[0], buildref.out, samples_ch, featureCountP)
      }
      else {
        fastqc(samples_ch)
        doSTAR(reference_path, samples_ch)
        FCounts(doSTAR.out[0],reference_path, samples_ch, featureCountP)
      }
      multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())

    emit:
      samples_ch
      qc_out = multiqc.out[1]
}

workflow {
  def date = new java.util.Date()
  def defaultReferenceId = "ensembl_v${params.ensembl_release ?: 107}_GRCh38"
  def reference_path = params.ref ?: "${params.reference_root ?: '/ref'}/${params.reference_id ?: defaultReferenceId}"
  params.reference_manifest = referenceManifestFromPath(reference_path)
  params.vep_cache = "${reference_path}/VEP"

  resource_source = channel.fromPath(params.scriptDir, checkIfExists: true)
  STAGE_WORKFLOW_RESOURCES(resource_source)

  // 1. Check_samples
  def sample_checked_csv
  if (params.check_samples != false) {
    def check_samples_res = Check_samples(
      params.csvSample,
      params.sampleInputDir,
      STAGE_WORKFLOW_RESOURCES.out.resources.map { resource_dir -> resource_dir.resolve('check_nb_sample.py') }
    )
    sample_checked_csv = check_samples_res[0]
    def sample_status = check_samples_res[1].map{ it.text.trim() }

    // 1.bis Vérifier samples
    sample_status
      .map { status ->
        if (status != "OK") {
          throw new RuntimeException("❌ Some sample fastqs are missing. Please check the sample list.")
        }
        return status
      }
      .set { verified_status_ch }
  }
  else {
    log.info("Skipping local sample directory scan. The CSV must include ID_Sample,suffix1,suffix2 columns.")
    sample_checked_csv = Channel.value(file(params.csvSample))
  }

  // 2. Construire le canal samples_ch depuis CSV
  def samples_ch

  if (params.single_end) {
    println "⚠️ Single-end mode activated"
    samples_ch = sample_checked_csv
    .splitCsv(header: true, sep:",")
    .map {
      [ it.ID_Sample,
        file("${params.sampleInputDir}/${it.ID_Sample}${it.suffix1}.fastq.gz")
      ]
    }
  } else {
    println "⚠️ Paired-end end mode activated"
    samples_ch = sample_checked_csv
    .splitCsv(header: true, sep:",")
    .map {
      [ it.ID_Sample,
        [
          file("${params.sampleInputDir}/${it.ID_Sample}${it.suffix1}.fastq.gz"),
          file("${params.sampleInputDir}/${it.ID_Sample}${it.suffix2}.fastq.gz")
        ]
      ]
    }
  }


  // 2bis. Extraire la liste des noms d’échantillons comme channel
  def list_names_ch = samples_ch.map { it[0] }.toList()


  // 3. Create_md5|Verify_md5 selon cas routine|reprocess
  def md5_error_file
  if (params.check_md5 != false) {
    if(params.routine==true){
      Create_md5(params.sampleInputDir, sample_checked_csv,
        STAGE_WORKFLOW_RESOURCES.out.resources.map { resource_dir -> resource_dir.resolve('createCheckSum.py') })
      md5_error_file = Create_md5.out[1]
    }else {
      Verify_md5(params.sampleInputDir, file(params.md5),
        STAGE_WORKFLOW_RESOURCES.out.resources.map { resource_dir -> resource_dir.resolve('verifyCheckSum.py') })
      md5_error_file = Verify_md5.out[1]
      // 3.bis Vérifier MD5 status file
      md5_error_file
        .map { it.text.trim() }
        .subscribe { status_md5 ->
        if (status_md5 != "OK") {
          throw new RuntimeException("❌ ERROR checking md5, please check your tab and fastq dir.")
          }
        }.set { verified_status_md5 }
    }
  }
  else {
    log.info("Skipping local MD5 scan for object-storage inputdir: ${params.sampleInputDir}")
    md5_error_file = Channel.value("SKIPPED")
  }

  // 4. Analyse RNAseq
  if (params.single_end) {
    Analysis_SE(samples_ch, md5_error_file, reference_path)
  }
  else {
    Analysis_PE(samples_ch, sample_checked_csv, params.sampleInputDir, md5_error_file, reference_path)
  }

  // 5. Vérifie les outputs process pour tous les échantillons.
  // Object-store output paths cannot be scanned with local filesystem APIs.
  def should_check_process = params.check_process_outputs != false && !params.outputdir.toString().startsWith('s3://')
  def process_status_ch
  if (should_check_process) {
    if (params.single_end) {
      process_status_ch = Check_process(sample_checked_csv, Analysis_SE.out.qc_out.collect(),
        STAGE_WORKFLOW_RESOURCES.out.resources.map { resource_dir -> resource_dir.resolve('check_process_nf.py') })[0]
    }
    else {
      process_status_ch = Check_process(sample_checked_csv, Analysis_PE.out.qc_out.collect(),
        STAGE_WORKFLOW_RESOURCES.out.resources.map { resource_dir -> resource_dir.resolve('check_process_nf.py') })[0]
    }
  }
  else {
    log.info("Skipping local output scan for object-storage outputdir: ${params.outputdir}")
    process_status_ch = Channel.value("OK")
  }

  process_status_ch
    .map { status_file -> should_check_process ? status_file.text.trim() : status_file }
    .subscribe { status_process ->
      if (status_process != "OK") {
        throw new RuntimeException("❌ Some process are missing, please check logs and outputs.")
      }
      log.info("✅ All process done successfully.")
    }

  list_names_ch
    .map { names ->
      [
        pipeline   : "RNApipeline NF - V1.7.0",
        pipeline_version: "1.7.0",
        date       : date.toString(),
        genome     : reference_path,
        reference_id: referenceIdFromPath(reference_path),
        reference_manifest: referenceManifestFromPath(reference_path),
        version    : makeVersionTag("V1.7", reference_path),
        sequenceID : params.runNumber,
        single_end : params.single_end,
        samples    : names,
        region_bed : params.bed,
        fastq_path : params.sampleInputDir,
        outputdir  : params.outputdir,
        scriptDir  : params.scriptDir,
        routine    : params.routine,
        run_star   : params.star,
        run_multiqc: params.multiqc,
        run_fastqc : params.fastqc,
        run_fcounts: params.fcounts,
        run_kallisto: params.kallisto,
        run_gatk   : params.gatk4,
        run_samtools : params.samtools_depth,
        run_deepvar: params.deepvariant,
        run_mpileup: params.mpileup,
        run_vep    : params.vep,
        containers : params.containers
      ]
    }
    .set { summary_ch }

  Write_summary(summary_ch)

}
