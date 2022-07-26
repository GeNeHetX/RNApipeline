
nextflow.enable.dsl=2


Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  [file(params.sampleInputDir + "/" + it + params.samPsuffix1) ,file(params.sampleInputDir + "/" + it + params.samPsuffix2)]] }
.set { samples_ch }


include {doSTAR; FCounts; multiqc} from '../modules/rna_seq_pipe.nf'
include {gatk_vc;Vep} from '../modules/variant_calling.nf'
include {Kallisto_paired_end} from '../modules/kallisto.nf'
include {buildref} from '../modules/index.nf'


workflow {
	
	if(params.ref=="no_ref") {
 	buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
	doSTAR(buildref.out, samples_ch)
	FCounts(doSTAR.out[0].collect(),buildref.out)
	Kallisto_paired_end(buildref.out, samples_ch)
	gatk_vc(doSTAR.out[0], buildref.out)
	}
	else {
	doSTAR(params.ref, samples_ch)
	FCounts(doSTAR.out[0].collect(),params.ref)
	Kallisto_paired_end(params.ref, samples_ch)
	gatk_vc(doSTAR.out[0], params.ref)
	}
	Vep(gatk_vc.out)
	multiqc(doSTAR.out[2].mix(doSTAR.out[1]).collect())
}
