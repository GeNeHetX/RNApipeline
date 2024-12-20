nextflow.enable.dsl=2


process Mosdepth {
    //publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1'

    input:
    // path("*${params.suffixbam}.bam")
    file(bamfile)
    file(baifile)
    tuple val(sample), file(fqFile)

	output:
        // path("${sample}_coverage")
    file("${sample}_coverage.per-base.bed.gz")
	// tuple val(sample), path('*.mosdepth.global.dist.txt')
    // tuple val(sample), path('*.mosdepth.summary.txt')
   // val(sample)
    //path('*.per-base.bed.gz')

    // tuple val(sample), path('*.per-base.bed.gz.csi')
    


	when:
	params.deepvariant == true
	"""
   
    mosdepth --threads ${params.nproc} \
    ${sample}_coverage  ${bamfile}
    

	"""
}


process Bedtools {
    //publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
        file(perbasebed)
        tuple val(sample), file(fqFile)
	//tuple val(sample), path('*.per-base.bed.gz')


	output:
    file("${sample}_3x.bed")


	when:
	params.deepvariant == true
	// """
    // gzip -dc ${perbasebed} | \
    // egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | \
    // awk -v OFS="\t" -v min_coverage=${params.min_coverage} '$4 >= min_coverage { print }' | \
    // bedtools merge -d 1 -c 4 -o mean -i - > ${sample}_3x.bed
	// """

    shell:
    '''
    gzip -dc !{perbasebed} | \
    egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | \
    awk -v OFS="\t" -v min_coverage=!{params.min_coverage} '$4 >= min_coverage { print }' | \
    bedtools merge -d 1 -c 4 -o mean -i - > !{sample}_3x.bed
	'''
}


process Deepvariant {
    publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'google/deepvariant:1.4.0'

	input :
    file(bamfile)
    file(baifile)
    tuple val(sample), file(fqFile)
    file(bedfile)
    path(reffa)
    path(model)

	output: 
    file("${sample}.output.vcf.gz")


	when:
	params.deepvariant == true
	"""
    run_deepvariant \
    --model_type=WES \
    --customized_model="${model}/model.ckpt" \
    --ref="${reffa}/ref.fa" \
    --reads="${bamfile}" \
    --output_vcf="${sample}.output.vcf.gz" \
    --num_shards=${params.nproc} \
    --regions=${bedfile} \
    --make_examples_extra_args="split_skip_reads=true,channels=''" \
    --intermediate_results_dir output/intermediate_results_dir
    
    """
    
}