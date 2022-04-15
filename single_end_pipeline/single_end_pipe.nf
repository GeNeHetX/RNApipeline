nextflow.enable.dsl=2




Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix1 )] }
.set { samples_ch}

include {doSTAR; FCounts; multiqc} from './rna_seq_pipe'

workflow {
	
	doSTAR(params.ref, samples_ch)
	FCounts(doSTAR.out[0],params.ref)
	multiqc(doSTAR.out[2],doSTAR.out[1])
	
}
