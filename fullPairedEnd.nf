
nextflow.enable.dsl=2

params.featureCountP=" -p "

Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  [file(params.sampleInputDir + "/" + it + params.samPsuffix1) ,file(params.sampleInputDir + "/" + it + params.samPsuffix2)]] }
  .set { samples_ch }


  include {doSTAR; samtools_index; FCounts; multiqc; fastqc; doOnlySTARnCount} from './modules/rna_seq_pipe.nf'
  include {gatk_vc;Vep} from './modules/variant_calling.nf'
  include {KallistoPE} from './modules/kallisto.nf'
  include {buildref} from './modules/index.nf'

log.info """\
RNAPIPE FULL PAIRED END - NF V1.6.0
===================================
genome : ${params.ref}
fastq : ${params.sampleInputDir}
results outputdir : ${params.outputdir}

run star: ${params.star}
run fcounts: ${params.fcounts}
run kallisto: ${params.kallisto}
run multiqc: ${params.multiqc}
run variant calling: ${params.variant_calling}
"""

  workflow {

    if(params.ref=="no_ref") {
        buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
        fastqc(samples_ch)
        doSTAR(buildref.out, samples_ch)
        FCounts(doSTAR.out[0].collect(),buildref.out)
        KallistoPE(buildref.out, samples_ch)
        gatk_vc(doSTAR.out[0], buildref.out)
      }
    else {
      fastqc(samples_ch)
      doSTAR(params.ref, samples_ch)
      samtools_index(doSTAR.out[0].collect())
      FCounts(doSTAR.out[0].collect(),params.ref)
      KallistoPE(params.ref, samples_ch)
      gatk_vc(doSTAR.out[0], params.ref)
    }
    
    // VCF annotation with VEP
    if (params.variant_calling == true){
      Vep(gatk_vc.out, params.ref)
    }

    //Agregate quality results
    multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
  }
