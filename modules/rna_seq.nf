nextflow.enable.dsl=2




Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix )] }
.set { samples_ch}

include {buildIndex; doSTAR; FCounts; multiqc} from './rna_seq_pipe'

workflow {
	buildIndex(params.refFasta, params.refGTF,params.vcf)
	doSTAR(buildIndex.out, samples_ch)
	FCounts(doSTAR.out[0],buildIndex.out)
	multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
	
}