nextflow.enable.dsl=2


params.featureCountP="  "

Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix1 )] }
.set { samples_ch}




  include {doOnlySTARnCount} from './modules/rna_seq_pipe.nf'

  include {buildref} from './modules/index.nf'


  workflow {
    if(params.ref== "no_ref" ){
      buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
      doOnlySTARnCount(buildref.out, samples_ch)

    }
    else {
      doOnlySTARnCount(params.ref, samples_ch)


    }
  }
