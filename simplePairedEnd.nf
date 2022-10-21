
nextflow.enable.dsl=2

params.featureCountP=" -p "

Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  [file(params.sampleInputDir + "/" + it + params.samPsuffix1) ,file(params.sampleInputDir + "/" + it + params.samPsuffix2)]] }
  .set { samples_ch }


  include {doSTAR; FCounts; multiqc; doOnlySTARnCount} from './modules/rna_seq_pipe.nf'
  include {Kallisto_paired_end} from './modules/kallisto.nf'
  include {buildref} from './modules/index.nf'


  workflow {


    if(params.ref== "no_ref" ){
      buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
      doOnlySTARnCount(buildref.out, samples_ch)
      Kallisto_paired_end(buildref.out, samples_ch)

    }
    else {
      doOnlySTARnCount(params.ref, samples_ch)
      Kallisto_paired_end(params.ref, samples_ch)

    }


  }
