nextflow.enable.dsl=2

params.ref_root = '/ref'
params.work_root = '/srv/slurm/scratch'
params.force = false
params.check_only = false
params.builder = "${projectDir}/ref_build_tod.sh"

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

process BUILD_REFERENCE {
    tag 'RNApipeline reference'
    label 'arm64'
    cpus 16
    memory '64 GB'
    time '48h'
    // ref_build_tod.sh manages its own node-local build directory.
    scratch false

    input:
    path builder, name: 'ref_build_tod.sh'

    output:
    path 'reference-build.done'

    script:
    def optional_args = [
        params.force ? '--force' : null,
        params.check_only ? '--check-only' : null
    ].findAll().join(' ')
    """
    bash ref_build_tod.sh \
        --ref-root '${params.ref_root}' \
        --work-root '${params.work_root}' \
        ${optional_args}
    touch reference-build.done
    """
}

workflow {
    builder_source = channel.fromPath(
        params.builder,
        checkIfExists: true
    )
    STAGE_REFERENCE_BUILDER(builder_source)
    BUILD_REFERENCE(STAGE_REFERENCE_BUILDER.out.builder)
}
