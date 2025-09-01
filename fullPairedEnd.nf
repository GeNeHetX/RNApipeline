
nextflow.enable.dsl=2

params.featureCountP=" -p "
//params.simpleCount = false


include {doSTAR; FCounts; multiqc; doOnlySTARnCount; samtools_index; fastqc} from './modules/rna_seq_pipe.nf'
include {gatk_vc; Vep as Vep_gatk; Vep as Vep_deepvariant; Vep as Vep_mpileup} from './modules/variant_calling.nf'
include {Mosdepth; Bedtools; Deepvariant} from './modules/deepvariant.nf'
include {bcftools_mpileup} from './modules/mpileup.nf'
include {KallistoPE} from './modules/kallisto.nf'
include {buildref} from './modules/index.nf'
include {Create_md5; Verify_md5; Check_samples; Check_process} from './modules/check_prepost_pipeline.nf'


workflow Analysis_PE{
    take: 
      samples_ch
      sample_csv
      fastqDir
    
    main:
        
      if(params.ref=="no_ref") {
        buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
        fastqc(samples_ch)
        doSTAR(buildref.out, samples_ch)
        FCounts(doSTAR.out[0], buildref.out, samples_ch)
        KallistoPE(buildref.out, samples_ch)
        gatk_vc(doSTAR.out.bam4bai, buildref.out)
      }
      else {
        fastqc(samples_ch)
        doSTAR(params.ref, samples_ch)
        samtools_index(doSTAR.out.bam4bai)
        KallistoPE(params.ref, samples_ch)

        if (params.fcounts == true){
          FCounts(doSTAR.out[0],params.ref, samples_ch)
        }

        // VC with gatk4 + vep
        if (params.gatk4 == true){
          gatk_vc(doSTAR.out.bam4bai, params.ref)
          Vep_gatk(gatk_vc.out.vc_file, params.ref,"gatk4")
        }

        // VC with mpileup + vep
        if (params.mpileup == true){
          bcftools_mpileup(samtools_index.out.align_files , params.ref, params.bed)
          Vep_mpileup(bcftools_mpileup.out.vc_file, params.ref,"mpileup")
        }

        // VC with deepvariant + vep
        if (params.deepvariant == true){
          Mosdepth(doSTAR.out[0],samtools_index.out[0],samples_ch)
          Bedtools(Mosdepth.out[0],samples_ch)
          Deepvariant(doSTAR.out[0], samtools_index.out[0], samples_ch, Bedtools.out[0], params.ref, params.modelckptdeepar)
          Vep_deepvariant(Deepvariant.out, params.ref,"deepvariant")
        }
       
      }
      //Agregate quality results
      if (params.multiqc == true){
        multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
      }

    emit:
      samples_ch
      qc_out = multiqc.out[1]

}

def date = new java.util.Date()
import groovy.json.JsonOutput

workflow Main {
  // 1. Check_samples
  def check_samples_res = Check_samples(params.csvSample, params.sampleInputDir, file("${params.scriptDir}/check_nb_sample.py"))
  def sample_checked_csv = check_samples_res[0]
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

  // 2. Construire le canal samples_ch depuis CSV
  def samples_ch = sample_checked_csv
    .splitCsv(header: true, sep:",")
    .map {
      [ it.ID_Sample,
        [
          file("${params.sampleInputDir}/${it.ID_Sample}${it.suffix1}.fastq.gz"),
          file("${params.sampleInputDir}/${it.ID_Sample}${it.suffix2}.fastq.gz")
        ]
      ]
    }
  // 2bis. Extraire la liste des noms d’échantillons comme channel
  def list_names_ch = samples_ch.map { it[0] }.toList()
  
    // 3. Create_md5|Verify_md5 selon cas routine|reprocess
    def md5_error_file
    if(params.routine==true){
      Create_md5( params.sampleInputDir, sample_checked_csv)
      md5_error_file = Create_md5.out[1]
    }else {
      Verify_md5(params.sampleInputDir, params.md5)
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

    // 4. Analyse RNAseq
    Analysis_PE(samples_ch, sample_checked_csv, params.sampleInputDir)

    // 5. Vérifie les outputs process pour tous les échantillons
    def check_process_res = Check_process(sample_checked_csv, Analysis_PE.out.qc_out.collect(),file("${params.scriptDir}/check_process_nf.py"))
    
    check_process_res[0]
      .map { it.text.trim() }
      .subscribe { status_process ->
        if (status_process == "OK") {
          log.info("✅ All process done successfully.")
          list_names_ch.subscribe { names ->
            def report = [
              pipeline   : "RNApipeline NF PAIRED END - V1.7.0",
              date       : date.toString(),
              genome     : params.ref,
              version    : "RNAv1.7.0_Ensemblv107",
              sequenceID : params.runNumber,
              samples    : names,
              region_bed : params.bed,
              fastq_path : params.sampleInputDir,
              outputdir  : params.outputdir,
              scriptDir  : params.scriptDir,
              routine    : params.routine,
              run_star   : params.star,
              run_fastqc : params.fastqc,
              run_fcounts: params.fcounts,
              run_kallisto: params.kallisto,
              run_multiqc: params.multiqc,
              run_gatk   : params.gatk4,
              run_deepvar: params.deepvariant,
              run_mpileup: params.mpileup,
              run_vep    : params.vep

            ]
            def outFile = new File("${params.outputdir}/${params.runNumber}_pipeline_summary.json")
            outFile.text = JsonOutput.prettyPrint(JsonOutput.toJson(report))
            log.info "✅ Pipeline summary écrit : ${outFile}"
          }
        } else {
          throw new RuntimeException("❌ Some process are missing, please check logs and outputs.")
        }
      }

}