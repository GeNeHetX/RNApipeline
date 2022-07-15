nextflow.enable.dsl=2




Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  file(params.sampleInputDir + "/" + it + params.samPsuffix1 )] }
.set { samples_ch}

include {doSTAR; FCounts; multiqc} from '../modules/rna_seq_pipe.nf'
include {kallisto_single_end} form '../modules/kallisto.nf'

workflow {
	doSTAR(params.ref, samples_ch)
	FCounts(doSTAR.out[0].collect(),params.ref)
	kallisto_single_end(params.ref, samples_ch)
	gatk_vc(doSTAR.out[0], params.ref)
	Vep(gatk_vc.out)
	multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
}
