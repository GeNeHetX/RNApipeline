nextflow.enable.dsl=2

params.ref_root = '/ref'
params.force = false
params.check_only = false
params.builder = "${projectDir}/ref_build_tod.sh"
params.kallisto_wrapper = "${projectDir}/../containers/kallisto/kallisto-wrapper.sh"
params.stage_root = null

process STAGE_REFERENCE_BUILDER {
    tag 'RNApipeline reference builder source'
    executor 'local'
    scratch false
    stageInMode 'copy'

    input:
    path builder_source, name: 'builder.source.sh'

    output:
    path 'ref_build_tod.sh', emit: builder

    script:
    """
    install -m 0755 builder.source.sh ref_build_tod.sh
    """
}

process PREPARE_REFERENCE_INPUTS {
    tag 'RNApipeline reference inputs'
    cpus 4
    memory '16 GB'
    time '12h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    val stage_root

    output:
    path 'reference-inputs.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase prepare-inputs \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch reference-inputs.ready
    """
}

process PREPARE_REFERENCE_IMAGES {
    tag 'RNApipeline reference containers'
    cpus 4
    memory '16 GB'
    time '4h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    val stage_root

    output:
    path 'reference-images.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase prepare-images \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch reference-images.ready
    """
}

process PUBLISH_KALLISTO_WRAPPER {
    tag 'RNApipeline Kallisto wrapper'
    executor 'local'
    scratch false

    input:
    path wrapper_source, name: 'kallisto-wrapper.sh'
    path images_ready
    val ref_root

    output:
    path 'kallisto-wrapper.ready', emit: ready

    script:
    """
    install -d -m 0775 '${ref_root}/tools/rnapipeline/v1.7.0/kallisto/0.51.1/bin'
    install -m 0755 kallisto-wrapper.sh \
        '${ref_root}/tools/rnapipeline/v1.7.0/kallisto/0.51.1/bin/kallisto'
    touch kallisto-wrapper.ready
    """
}

process BUILD_SEQUENCE_INDEXES {
    tag 'RNApipeline sequence indexes'
    cpus 16
    memory '32 GB'
    time '8h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    val stage_root

    output:
    path 'sequence-indexes.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase sequence-indexes \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch sequence-indexes.ready
    """
}

process BUILD_VARIANT_INDEX {
    tag 'RNApipeline variant index'
    cpus 8
    memory '16 GB'
    time '4h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    val stage_root

    output:
    path 'variant-index.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase variant-index \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch variant-index.ready
    """
}

process BUILD_KALLISTO_INDEX {
    tag 'RNApipeline Kallisto index'
    cpus 16
    memory '32 GB'
    time '8h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    val stage_root

    output:
    path 'kallisto-index.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase kallisto-index \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch kallisto-index.ready
    """
}

process BUILD_GTF_ARTIFACTS {
    tag 'RNApipeline GTF artifacts'
    cpus 8
    memory '16 GB'
    time '4h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    val stage_root

    output:
    path 'gtf-artifacts.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase gtf-artifacts \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch gtf-artifacts.ready
    """
}

process BUILD_STAR_INDEX {
    tag 'RNApipeline STAR index'
    cpus 16
    memory '64 GB'
    time '24h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    val stage_root

    output:
    path 'star-index.ready', emit: ready

    script:
    """
    bash ref_build_tod.sh \\
        --phase star-index \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch star-index.ready
    """
}

process CHECK_REFERENCE {
    tag 'RNApipeline reference preflight'
    cpus 1
    memory '2 GB'
    time '30m'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    val stage_root

    output:
    path 'reference-preflight.ready'

    script:
    """
    bash ref_build_tod.sh \\
        --check-only \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}'
    touch reference-preflight.ready
    """
}

process FINALIZE_REFERENCE {
    tag 'RNApipeline reference publication'
    cpus 16
    memory '64 GB'
    time '12h'
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'
    path inputs_ready
    path images_ready
    path sequence_ready
    path variant_ready
    path kallisto_ready
    path gtf_ready
    path star_ready
    val stage_root

    output:
    path 'reference-build.done'

    script:
    def force_arg = params.force ? '--force' : ''
    """
    bash ref_build_tod.sh \\
        --phase finalize \\
        --ref-root '${params.ref_root}' \\
        --work-root '${stage_root}' \\
        ${force_arg}
    touch reference-build.done
    """
}

workflow {
    builder_source = channel.fromPath(params.builder, checkIfExists: true)
    staged = STAGE_REFERENCE_BUILDER(builder_source)
    stage_root = params.stage_root ?: "${workflow.workDir}/reference-stage"

    if (params.check_only) {
        CHECK_REFERENCE(staged, stage_root)
    } else {
        inputs = PREPARE_REFERENCE_INPUTS(staged, stage_root)
        images = PREPARE_REFERENCE_IMAGES(staged, stage_root)
        wrapper_source = channel.fromPath(
            params.kallisto_wrapper, checkIfExists: true
        )
        PUBLISH_KALLISTO_WRAPPER(
            wrapper_source, images, params.ref_root
        )

        sequence = BUILD_SEQUENCE_INDEXES(
            staged, inputs, images, stage_root
        )
        variants = BUILD_VARIANT_INDEX(
            staged, inputs, images, stage_root
        )
        kallisto = BUILD_KALLISTO_INDEX(
            staged, inputs, images, stage_root
        )
        gtf = BUILD_GTF_ARTIFACTS(
            staged, inputs, images, stage_root
        )
        star = BUILD_STAR_INDEX(
            staged, inputs, images, stage_root
        )

        FINALIZE_REFERENCE(
            staged,
            inputs,
            images,
            sequence,
            variants,
            kallisto,
            gtf,
            star,
            stage_root
        )
    }
}
