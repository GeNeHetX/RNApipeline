nextflow.enable.dsl=2




Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix1 )] }
.set { samples_ch}

include doSTAR; FCounts; multiqc} from '../modules/rna_seq_pipe.nf'

workflow {
	doSTAR(buildIndex.out, samples_ch)
	FCounts(doSTAR.out[0],params.ref)
	multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
}
