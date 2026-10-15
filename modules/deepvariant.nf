nextflow.enable.dsl=2


process Mosdepth {
    label 'amd64'
    //publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1'

    input:
    // path("*${params.suffixbam}.bam")
    tuple val(sample), path(bamfile), path(baifile)

	output:
        // path("${sample}_coverage")
    tuple val(sample), path("${sample}_coverage.per-base.bed.gz"), emit: coverage
	// tuple val(sample), path('*.mosdepth.global.dist.txt')
    // tuple val(sample), path('*.mosdepth.summary.txt')
   // val(sample)
    //path('*.per-base.bed.gz')

    // tuple val(sample), path('*.per-base.bed.gz.csi')



	when:
	params.deepvariant.toString().toBoolean()
	script:
	"""

    mosdepth --threads ${task.cpus} \
    ${sample}_coverage  ${bamfile}


	"""
}


process Bedtools {
    label 'amd64'
    //publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
        tuple val(sample), path(perbasebed)
	//tuple val(sample), path('*.per-base.bed.gz')


	output:
    tuple val(sample), path("${sample}_3x.bed"), emit: regions


	when:
	params.deepvariant.toString().toBoolean()
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
    label 'amd64'
    publishDir "${params.outputdir}/deepvariant_output", mode: 'copy'

    container 'google/deepvariant:1.4.0'

	input :
    tuple val(sample), path(bamfile), path(baifile), path(bedfile)
    path(reffa)
    path(model)

	output:
    tuple val(sample), file("${sample}.output.vcf.gz")


	when:
	params.deepvariant.toString().toBoolean()
	script:
	"""
    run_deepvariant \
    --model_type=WES \
    --customized_model="${model}/model.ckpt" \
    --ref="${reffa}/ref.fa" \
    --reads="${bamfile}" \
    --output_vcf="${sample}.output.vcf.gz" \
    --num_shards=${task.cpus} \
    --regions=${bedfile} \
    --make_examples_extra_args="split_skip_reads=true,channels=''" \
    --intermediate_results_dir output/intermediate_results_dir

    """

}
