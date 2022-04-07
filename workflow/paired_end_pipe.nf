
nextflow.enable.dsl=2


Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  [file(params.sampleInputDir + "/" + it + params.samPsuffix1) ,file(params.sampleInputDir + "/" + it + params.samPsuffix2)]] }
.set { samples_ch }


include {buildIndex; doSTAR; FCounts; multiqc} from '/home/nassimaima/Rna-seq_pipeline-main/modules/rna_seq_pipe.nf'
include {PicardSamsorted; gatk_vc;Vep} from '/home/nassimaima/Rna-seq_pipeline-main/modules/variant_calling.nf'

workflow {
	buildIndex(params.refFasta, params.refGTF,params.vcf)
	doSTAR(buildIndex.out, samples_ch)
	FCounts(doSTAR.out[0],buildIndex.out)
	PicardSamsorted(doSTAR.out[0])
	gatk_vc(buildIndex.out,PicardSamsorted.out)
	Vep(gatk_vc.out)
	multiqc(doSTAR.out[2],doSTAR.out[1])
}
