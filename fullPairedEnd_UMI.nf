
nextflow.enable.dsl=2

params.featureCountP=" -p "

Channel.fromList(file(params.sampleList).readLines())
.map {
  [it ,  [file(params.sampleInputDir + "/" + it + params.samPsuffix1) ,file(params.sampleInputDir + "/" + it + params.samPsuffix2)]] }
  .set { samples_ch }


  include {doSTAR; samtools_index; samtools_depth; FCounts; multiqc; fastqc; doOnlySTARnCount} from './modules/rna_seq_pipe.nf'
  include {gatk_vc; Vep as Vep_gatk; Vep as Vep_deepvariant; Vep as Vep_mpileup} from './modules/variant_calling.nf'
  include {KallistoPE} from './modules/kallisto.nf'
  include {remove_multimapped_reads; bcftools_mpileup} from './modules/mpileup.nf'
  include {buildref} from './modules/index.nf'
  include {Mosdepth; Bedtools; Deepvariant} from './modules/deepvariant.nf'
  include {UMItools_extract} from './modules/UMI.nf' //UMItools_dedup

log.info """\
RNAPIPE FULL PAIRED END - NF V1.6.3
===================================
genome : ${params.ref}
fastq path : ${params.sampleInputDir}
samples list : ${params.sampleList}
region bed : ${params.bed}
results outputdir : ${params.outputdir}

run star: ${params.star}
run fastqc: ${params.star}
run fcounts: ${params.fcounts}
run samtools depth: ${params.samtools_depth}
run kallisto: ${params.kallisto}
run gatk4: ${params.gatk4}
run deepvariant: ${params.deepvariant}
run mpileup: ${params.mpileup}
remove multimapped: ${params.no_multimapped}
run vep: ${params.vep}
run multiqc: ${params.multiqc}
"""

workflow RNApipe_with_UMI {
  
}

  workflow {

    if(params.ref=="no_ref") {
        buildref(params.fasta_ref,params.GTF,params.cdna,params.known_vcf)
        fastqc(samples_ch)
        doSTAR(buildref.out, samples_ch)
        FCounts(doSTAR.out[0],buildref.out, samples_ch)
        KallistoPE(buildref.out, samples_ch)
        gatk_vc(doSTAR.out[0], buildref.out)
      }
    else {
      fastqc(samples_ch)
      doSTAR(params.ref, samples_ch)
      samtools_index(doSTAR.out.bam4bai)
      KallistoPE(params.ref, samples_ch)
      if (params.fcounts == true){
        FCounts(doSTAR.out[0],params.ref, samples_ch)
      }

      if (params.samtools_depth == true){
        samtools_depth(doSTAR.out.bam4bai, params.bed)
      }

      if (params.UMI == true){
        UMItools_extract(samples_ch)
        UMItools_dedup(doSTAR.out.bam4bai)
      }

      // VC with gatk4 + vep
      if (params.gatk4 == true){
        gatk_vc(doSTAR.out.bam4bai, params.ref)
        Vep_gatk(gatk_vc.out.vc_file, params.ref,"gatk4")
      }

      // VC with mpileup + vep
      if (params.mpileup == true & params.no_multimapped == true){
        remove_multimapped_reads(samtools_index.out.align_files)
        bcftools_mpileup(remove_multimapped_reads.out.filtered_bam , params.ref, params.bed)
        Vep_mpileup(bcftools_mpileup.out.vc_file, params.ref,"mpileup")
      }

      else if (params.mpileup == true){ // no_multimapped == false
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
  }
